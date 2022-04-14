#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 22:26:00 2021

@author: firat
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymongo import MongoClient

surfen_crystalium = {'mp-149': {'322': 1.39,
                                '320': 1.45,
                                '210': 1.49,
                                '310': 1.49,
                                '321': 1.52,
                                '311': 1.52,
                                '111': 1.57,
                                '331': 1.58,
                                '332': 1.61,
                                '221': 1.65,
                                '211': 1.67,
                                '110': 1.76,
                                '100': 2.13},
                     'mp-134': {'111': 0.77,
                                '322': 0.90,
                                '100': 0.91,
                                '332': 0.91,
                                '221': 0.95,
                                '331': 0.96,
                                '321': 0.96,
                                '311': 0.98,
                                '110': 0.98,
                                '211': 0.98,
                                '310': 0.99,
                                '320': 1.01,
                                '210': 1.02},
                     'mp-30': {'111': 1.34,
                               '332': 1.43,
                               '322': 1.45,
                               '100': 1.47,
                               '221': 1.48,
                               '331': 1.52,
                               '311': 1.54,
                               '110': 1.56,
                               '321': 1.58,
                               '310': 1.59,
                               '210': 1.60,
                               '320': 1.62,
                               '211': 1.63},
                     'mp-13': {'110': 2.45,
                               '100': 2.50,
                               '320': 2.56,
                               '210': 2.57,
                               '211': 2.61,
                               '311': 2.63,
                               '331': 2.63,
                               '321': 2.63,
                               '221': 2.66,
                               '332': 2.68,
                               '322': 2.68,
                               '111': 2.73,
                               '310': 3.43}}


# with MPRester() as m:
#     struct = m.get_structure_by_material_id('mp-13', conventional_unit_cell=True)
# millers = get_symmetrically_distinct_miller_indices(struct, 3)
# cryst_data = {}
# with MPRester() as m:
#     for miller in millers:
#         surfen = m.get_surface_data('mp-13', miller)['surface_energy']
#     millerstr = ''.join([str(i) for i in miller])
#     cryst_data[millerstr] = surfen

def get_surfen_dict_from_mpid(mpid, coll, fig_str=None, db_file='auto', high_level=True, to_poscar=False):
    client = MongoClient(
        "mongodb+srv://tfroot:letsseeifthispasswordissecureenough@triboflow.pzgj2.mongodb.net/myFirstDatabase?retryWrites=true&w=majority")
    tf_surfen_coll = client['surfen_test']['PBE.slab_data.LEO']
    data = tf_surfen_coll.find_one({'mpid': f'{mpid}'})

    if not fig_str:
        fig_str = f'{mpid}_{coll}_surfen.png'

    tf_bulk_coll = client['surfen_test']['PBE.bulk_data']
    bulk_data = tf_bulk_coll.find_one({'mpid': mpid})
    bulk_prim = Structure.from_dict(bulk_data['structure_equiVol'])
    bulk_conv = sga(bulk_prim).get_conventional_standard_structure()
    lattice = bulk_conv.lattice

    surfen_dict = {}
    surfen_data = data['miller_list']
    for hkl, uids in surfen_data.items():
        hkl_list = [int(a) for a in hkl]
        surfen_dict[hkl] = {}
        surfen_by_hkl = []
        for uid, data in uids.items():
            try:
                surfen = data['surface_energy']
            except KeyError:
                continue
            uid_short = uid[:2]
            min_surfen = min(surfen['top'], surfen['bottom'])
            surfen_by_hkl.append(min_surfen)
            if to_poscar:
                slab = Slab.from_dict(data['structure'])
                formula = slab.composition.reduced_formula
                slab.to('poscar', f'./poscars/{formula}-{hkl}_{uid_short}-TF.vasp')
        min_surfen_hkl = min(surfen_by_hkl)
        surfen_dict[hkl] = min_surfen_hkl

    color_n = len(surfen_dict)
    colors = sns.color_palette("Spectral", color_n).as_hex()

    fig, ax = plt.subplots()
    fig.dpi = 200
    count = 0
    for k, v in surfen_dict.items():
        if k == "310" and mpid == 'mp-13':
            continue
        surfen_lit = surfen_crystalium[mpid][k]
        surfen_our = v
        ax.scatter(surfen_lit, surfen_our, c=colors[count], label=k, s=100, edgecolors='black')
        ax.text(surfen_lit, surfen_our, k)
        count += 1
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel(r'$\gamma^{crystalium}$ [$J/m^2$]', fontsize=20)
    ax.set_ylabel(r'$\gamma^{tw}$ [$J/m^2$]', fontsize=20)
    ax.locator_params(tight=True, nbins=4)
    ax.legend(bbox_to_anchor=(1.025, 1), prop={'size': 8})
    # fig.savefig(fig_str, bbox_inches='tight', dpi=300)

    return lattice, surfen_dict

title = st.text_input('Materials Project ID', 'mp-134')


try:
    lattice, surfen_dict = get_surfen_dict_from_mpid(title,
                                                 coll='PBE.slab_data.LEO',
                                                 high_level='surfen_test',
                                                 to_poscar=False)
except:
    st.write('No results for', title)
else:
    miller_list = [[int(x) for x in a] for a in list(surfen_dict.keys())]
    e_surf_list = list(surfen_dict.values())
    wulff = WulffShape(lattice, miller_list, e_surf_list)
    wulff_plotly = wulff.get_plotly()
    # wulff_plotly.show()

    surfen_dict = pd.DataFrame(surfen_dict.values(), surfen_dict.keys())

    st.write('Showing results for', title)
    st.dataframe(surfen_dict)

    st.plotly_chart(wulff_plotly)
