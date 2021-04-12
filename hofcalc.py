import streamlit as st
import pubchempy as pcp
import math
import pandas as pd
from pyvalem.formula import Formula
from collections import defaultdict

volumes = {"H" : 5.08, "He" : math.nan, "Li" : 22.6, "Be" : 36., "B" : 13.24,
    "C" : 13.87, "N" : 11.8, "O" : 11.39, "F" : 11.17, "Ne" : math.nan,
    "Na" : 26.0, "Mg" : 36., "Al" : 39.6, "Si" : 37.3, "P" : 29.5, "S" : 25.2,
    "Cl" : 25.8, "Ar" : math.nan, "K" : 36., "Ca" : 45., "Sc" : 42.,
    "Ti" : 27.3, "V" : 24., "Cr" : 28.1, "Mn" : 31.9, "Fe" : 30.4, "Co" : 29.4,
    "Ni" : 26., "Cu" : 26.9, "Zn" : 39., "Ga" : 37.8, "Ge" : 41.6, "As" : 36.4,
    "Se" : 30.3, "Br" : 32.7, "Kr" : math.nan, "Rb" : 42., "Sr" : 47.,
    "Y" : 44., "Zr" : 27., "Nb" : 37., "Mo" : 38., "Tc" : 38, "Ru" : 37.3,
    "Rh" : 31.2, "Pd" : 35., "Ag" : 35., "Cd" : 51., "In" : 55, "Sn" : 52.8,
    "Sb" : 48., "Te" : 46.7, "I" : 46.2, "Xe" : 45., "Cs" : 46., "Ba" : 66,
    "La" : 58.,    "Ce" : 54., "Pr" : 57., "Nd" : 50., "Pm" : math.nan,
    "Sm" : 50., "Eu" : 53., "Gd" : 56., "Tb" : 45., "Dy" : 50., "Ho" : 42.,
    "Er" : 54., "Tm" : 49., "Yb" : 59., "Lu" : 35., "Hf" : 40., "Ta" : 43.,
    "W" : 38.8, "Re" : 42.7, "Os" : 41.9, "Ir" : 34.3, "Pt" : 38., "Au" : 43.,
    "Hg" : 38., "Tl" : 54., "Pb" : 52., "Bi" : 60., "Po" : math.nan,
    "At" : math.nan, "Rn" : math.nan, "Fr" : math.nan, "Ra" : math.nan,
    "Ac" : 74., "Th" : 56., "Pa" : 60., "U" : 58., "Np" : 45., "Pu" : math.nan,
    "Am" : 17., "Cm" : math.nan, "Bk" : math.nan, "Cf" : math.nan,
    "Es" : math.nan, "Fm" : math.nan}

def get_volume(formula, volumes=volumes, temperature=298):
    volume = 0
    alpha = 0.95e-4
    eighteen_angstrom_volume = 0
    for element in formula.keys():
        if element not in volumes:
            st.write("Element",element,"not recognized.")
            return None, None
        else:
            if element != "H":
                eighteen_angstrom_volume += 18*formula[element]
            if volumes[element] is not math.nan:
                volume += (volumes[element] *
                                formula[element])*(1+alpha*(temperature - 298))
            else:
                st.write("Element", element, "does not have a Hofmann volume")
                st.write("Using the 18-angstrom rule to estimate missing volume")
                volume += 18 * formula[element]
    return round(volume, 1), eighteen_angstrom_volume

def get_formula(autodetect=True, key=1):
    #formula_option = st.radio("",
    #        ["Search PubChem for compound","Enter chemical formula"], key=key)
    #if formula_option == "Enter chemical formula":
    col1, col2 = st.beta_columns([3,1])
    with col2:
        search_mode = st.selectbox("Search by:",
                    ["Name","SMILES", "InChI"], key=key)
    with col1:
        if search_mode == "Name":
            user_input = st.text_input("Enter chemical formula or compound name",
                    value='', max_chars=None,
                    key=key, type='default')
        elif search_mode == "SMILES":
            user_input = st.text_input("Enter chemical formula or compound SMILES",
                    value='', max_chars=None,
                    key=key, type='default')
        else:
            user_input = st.text_input("Enter chemical formula or compound InChI",
                    value='', max_chars=None,
                    key=key, type='default')
    if user_input != "":
        molecular_formulae = []
        if autodetect:
            components = user_input.split(" ")
        else:
            components = [user_input]
        for component in components:
            if component != "":
                try:
                    component = component.replace("|", "")
                    f = Formula(component)
                    mf = f.atom_stoich
                except:
                    try:
                        molecule = pcp.get_compounds(component, search_mode.lower())[0]
                        f = Formula(molecule.molecular_formula)
                        mf = f.atom_stoich
                    except:
                        st.write("Error")
                        st.write("Unable to parse",component,"as a chemical formula \
                        or find it on PubChem when searching by",search_mode+".")
                        return None
                molecular_formulae.append(mf)
        molecular_formula = defaultdict(int)
        for mf in molecular_formulae:
            for key, value in mf.items():
                molecular_formula[key] += value
    else:
        molecular_formula = None

    return molecular_formula


st.set_page_config(page_title='HofCalc WebApp', page_icon = None,
                    layout = 'centered', initial_sidebar_state = 'auto')


# Top section
st.title("HofCalc")
st.markdown("*Crystallographic Volume Estimation*")
st.sidebar.title("Options")

with st.sidebar:
    option = st.radio("Select function",
                ["Volume estimation","Display Hofmann volumes"])
    for i in range(6):
        st.write("")
    with st.beta_expander(label="References", expanded=False):
        st.write("Hofmann, D.W.M. (2002), Fast estimation of crystal densities.\
            Acta Cryst. B, 58: 489-493. \
            https://doi.org/10.1107/S0108768101021814")
        st.sidebar.write("")
        st.sidebar.write("")
        st.write("WebApp designed by Mark Spillman")

if option == "Volume estimation":
    mfs = []
    col1, col2, col3 = st.beta_columns(3)

    with col1:
        temperature = st.number_input("Temperature / K",
                                min_value=0, max_value=None, value=298, step=10)
    with col2:
        unit_cell_volume = st.number_input("Unit cell volume (optional)",
                                min_value=0, max_value=None, value=0, step=100)
    autodetect = st.checkbox("Autodetect fragments", value=True,
                                                                    key=None)
    with col3:
        if not autodetect:
            number_of_fragments = st.number_input("Number of fragments", min_value=1,
                                max_value=None, value=1, step=1)
        else:
            number_of_fragments = 1
    for i in range(number_of_fragments):
        if number_of_fragments > 1:
            st.markdown("***Fragment "+str(i+1)+"***")
        molecular_formula = get_formula(autodetect=autodetect, key=i)
        if molecular_formula is not None:
            mfs.append(molecular_formula)
    total_hofmann = 0
    total_18 = 0
    st.write("")
    st.write("")
    if len(mfs) != 0:
        st.markdown("**Results**")
        if number_of_fragments > 1:
            col1, col2 = st.beta_columns(2)
            with col1:
                expander1 = st.beta_expander(label="Individual Hofmann",
                                            expanded=False)
            with col2:
                expander2 = st.beta_expander(label="Individual 18 ångström",
                                            expanded=False)
        for i, molecular_formula in enumerate(mfs):
            if molecular_formula is not None:
                volume = get_volume(molecular_formula, temperature=temperature)
                if number_of_fragments > 1:
                    with expander1:
                        st.markdown("Hofmann volume (frag "+str(i+1)+
                                    ") = " + str(volume[0])+" $Å^3$")
                    with expander2:
                        st.markdown("18 ångström rule (frag "+str(i+1)+
                                    ") = " + str(volume[1])+" $Å^3$")
                total_hofmann += volume[0]
                total_18 += volume[1]
        if number_of_fragments > 1:
            st.write("")
            st.write("")
        col1, col2 = st.beta_columns(2)
        with col1:
            st.markdown("Total Hofmann volume = "+str(total_hofmann)+" $Å^3$")
        with col2:
            st.markdown("Total 18 ångström rule = "+str(total_18)+" $Å^3$")
        if unit_cell_volume != 0:
            with col1:
                st.write("Unit cell ÷ Hofmann =", str(round(unit_cell_volume / total_hofmann, 2)))
            with col2:
                st.write("Unit cell ÷ 18Å =", str(round(unit_cell_volume / total_18, 2)))

else:
    st.write("Average crystallographic volumes reported by Hofmann at 298 K \
            (see references)")
    volume_df = pd.DataFrame.from_dict(volumes, columns=["Volume"],
                                        orient="index", dtype=float)
    st.table(volume_df.style.format('{:.2f}'))


