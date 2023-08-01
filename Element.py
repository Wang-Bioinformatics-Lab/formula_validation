#!/usr/bin/env python

"""
This Python module contains not only the class Element, but also the test of
this Element class.

@contents :  This Python module contains not only the class Element, but also the test of
this Element class.
@project :  N/A
@program :  N/A
@file :  Element.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)

@version :  0.0.1, 20 July 2023
@information :  The Zen of Python
                  https://www.python.org/dev/peps/pep-0020/
                Style Guide for Python Code
                  https://www.python.org/dev/peps/pep-0008/
                Example NumPy Style Python Docstrings
                  http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
                doctest – Testing through documentation
                  https://pymotw.com/2/doctest/

@copyright :  Copyright 2023 GNU AFFERO GENERAL PUBLIC.
              All rights are reserved. Reproduction in whole or in part is
              prohibited without the written consent of the copyright owner.
"""
__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"

from enum import Enum

Element_type = Enum('Element_type',['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Ni', 'Co', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Nh', 'Fl', 'Uup', 'Mc', 'Lv', 'Uus', 'Ts', 'Uuo'])
element_weights={Element_type.H:1.0078250321, Element_type.He:2.014102, Element_type.He:4.0026032542, Element_type.Li:7.016004558, Element_type.Be:9.012182, Element_type.B:11.009305, Element_type.C:12, Element_type.N:14.003074, Element_type.O:15.994915, Element_type.F:18.998403163, Element_type.Ne:19.99244, Element_type.Na:22.98976928, Element_type.Mg:23.98504, Element_type.Al:26.9815385, Element_type.Si:27.97693, Element_type.P:30.973761632, Element_type.S:31.9720710015, Element_type.Cl:34.96885, Element_type.Ar:39.962383, Element_type.K:38.96371, Element_type.Ca:39.96259, Element_type.Sc:44.95591, Element_type.Ti:47.94794, Element_type.V:50.94396, Element_type.Cr:51.94051, Element_type.Mn:54.93804, Element_type.Fe:55.93494, Element_type.Co:58.93319, Element_type.Ni:57.93534, Element_type.Cu:62.92960, Element_type.Zn:63.92914, Element_type.Ga:68.92557, Element_type.Ge:73.92118, Element_type.As:74.92159, Element_type.Se:79.91652, Element_type.Br:78.91834, Element_type.Kr:83.91150, Element_type.Rb:84.91179, Element_type.Sr:87.90561, Element_type.Y:88.90584, Element_type.Zr:89.90470, Element_type.Nb:92.90637, Element_type.Mo:97.90540, Element_type.Tc:98.0, Element_type.Ru:101.90434, Element_type.Rh:102.90550, Element_type.Pd:105.90348, Element_type.Ag:106.90509, Element_type.Cd:113.90337, Element_type.In:114.90388, Element_type.Sn:119.90220, Element_type.Sb:120.90381, Element_type.Te:129.90622, Element_type.I:126.90447, Element_type.Xe:131.90416, Element_type.Cs:132.90545, Element_type.Ba:137.90525, Element_type.La:138.90636, Element_type.Ce:139.90544, Element_type.Pr:140.90766, Element_type.Nd:141.90773, Element_type.Pm:145.0, Element_type.Sm:151.91974, Element_type.Eu:152.92124, Element_type.Gd:157.92411, Element_type.Tb:158.92535, Element_type.Dy:163.92918, Element_type.Ho:164.93033, Element_type.Er:165.93030, Element_type.Tm:168.93422, Element_type.Yb:173.93887, Element_type.Lu:174.94078, Element_type.Hf:179.94656, Element_type.Ta:180.94800, Element_type.W:183.95093, Element_type.Re:186.95575, Element_type.Os:191.96148, Element_type.Ir:192.96292, Element_type.Pt:194.96479, Element_type.Au:196.96657, Element_type.Hg:201.97064, Element_type.Tl:204.97443, Element_type.Pb:207.97665, Element_type.Bi:208.98040, Element_type.Po:209.98287, Element_type.At:210.0, Element_type.Rn:222.0, Element_type.Fr:223.0, Element_type.Ra:226.0, Element_type.Ac:227.0, Element_type.Th:232.03806, Element_type.Pa:231.03588, Element_type.U:238.05079, Element_type.Np:237.0, Element_type.Pu:244.0, Element_type.Am:243.0, Element_type.Cm:247.0, Element_type.Bk:247.0, Element_type.Cf:251.0, Element_type.Es:252.0, Element_type.Fm:257.0, Element_type.Md:258.0, Element_type.No:259.0, Element_type.Lr:266.0, Element_type.Rf:267.0, Element_type.Db:268.0, Element_type.Sg:269.0, Element_type.Bh:270.0, Element_type.Hs:269.0, Element_type.Mt:278.0, Element_type.Ds:281.0, Element_type.Rg:282.0, Element_type.Cn:285.0, Element_type.Uut:286.0, Element_type.Fl:289.0, Element_type.Uup:289.0, Element_type.Lv:293.0, Element_type.Uus:294.0, Element_type.Uuo:294.0}


def main():
    pass
    

if __name__ == "__main__":
    main()