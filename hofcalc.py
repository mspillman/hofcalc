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
    volume = 0.0
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
                                formula[element])*(1.+alpha*(temperature - 298.))
            else:
                st.write("Element", element, "does not have a Hofmann volume")
                st.write("Using the 18-angstrom rule to estimate missing volume")
                volume += 18 * formula[element]
    return round(volume, 2), eighteen_angstrom_volume

def get_formula(autodetect=True, key=1):
    #formula_option = st.radio("",
    #        ["Search PubChem for compound","Enter chemical formula"], key=key)
    #if formula_option == "Enter chemical formula":
    col1, col2 = st.beta_columns([3,1])
    with col2:
        search_mode = st.selectbox("Search by:",
                    ["Name","SMILES"], key=key)
    with col1:
        if search_mode == "Name":
            user_input = st.text_input("Enter chemical formula or compound name",
                    value='', max_chars=None,
                    key=key, type='default')
        elif search_mode == "SMILES":
            user_input = st.text_input("Enter chemical formula or compound SMILES",
                    value='', max_chars=None,
                    key=key, type='default')
    if user_input != "":
        if "<" in user_input:
            user_input = user_input.split("</font>")[0]
            user_input = user_input.split("<b>")[1]
            user_input = user_input.split(">")[1:]
            user_input = ("".join(
                                user_input)).replace("<sub","").replace(
                                                    "</sub","").replace(" ","")
            st.write("Interpreting formula as", user_input)
        molecular_formulae = []
        user_input = user_input.replace("+", "")
        user_input = user_input.replace("-", "")
        if autodetect:
            if "," in user_input:
                components = user_input.split(",")
                components = [x.replace(" ", "") for x in components]
            else:
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
                        other_modes = ["name", "smiles"]
                        other_modes.remove(search_mode.lower())
                        result = False
                        for om in other_modes:
                            try:
                                molecule = pcp.get_compounds(component, om)[0]
                                f = Formula(molecule.molecular_formula)
                                mf = f.atom_stoich
                                result = True
                                break
                            except:
                                pass
                        if not result:
                            st.write("Error")
                            st.write("Unable to parse",component,"as a chemical \
                                formula or find it on PubChem when searching by",
                                search_mode+".")
                            st.write("Subsequent search by",other_modes[0],"also failed")
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
                ["Volume estimation","Examples","Display Hofmann volumes"])
    for i in range(6):
        st.write("")
    with st.beta_expander(label="References", expanded=False):
        st.write("Hofmann, D.W.M. (2002), Fast estimation of crystal densities.\
            Acta Cryst. B, 58: 489-493. \
            https://doi.org/10.1107/S0108768101021814")
        st.write("")
        st.write("")
        st.write("PubChem: https://pubchem.ncbi.nlm.nih.gov/")
        st.write("")
        st.write("")
        st.write("WebApp designed by Mark Spillman")
        st.write("https://github.com/mspillman/hofcalc")

if option == "Volume estimation":
    mfs = []
    col1, col2, col3 = st.beta_columns(3)

    with col1:
        temperature = st.number_input("Temperature / K",
                                min_value=0, max_value=None, value=298, step=10)
    with col2:
        unit_cell_volume = st.number_input("Unit cell volume (optional)",
                                min_value=0.0, max_value=None, value=0.0,
                                step=100.0)
    autodetect = st.checkbox("Autodetect fragments", value=True,
                                                                    key=None)
    with col3:
        if not autodetect:
            number_of_fragments = st.number_input("Number of fragments",
                                min_value=1, max_value=None, value=1, step=1)
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
            #col1, col2 = st.beta_columns(2)
            #with col1:
            expander1 = st.beta_expander(label="Individual Hofmann",
                            expanded=False)
            #with col2:
            #    expander2 = st.beta_expander(label="Individual 18 ångström",
            #                                expanded=False)
        for i, molecular_formula in enumerate(mfs):
            if molecular_formula is not None:
                volume = get_volume(molecular_formula, temperature=temperature)
                if number_of_fragments > 1:
                    with expander1:
                        st.markdown("Hofmann volume (frag "+str(i+1)+
                                    ") = " + str(volume[0])+" $Å^3$")
                    #with expander2:
                    #    st.markdown("18 ångström rule (frag "+str(i+1)+
                    #                ") = " + str(volume[1])+" $Å^3$")
                total_hofmann += volume[0]
                total_18 += volume[1]
        if number_of_fragments > 1:
            st.write("")
            st.write("")
        col1, col2 = st.beta_columns(2)
        with col1:
            st.markdown("Total Hofmann volume = "+str(round(total_hofmann, 3))+" $Å^3$")
        #with col2:
        #    st.markdown("Total 18 ångström rule = "+str(total_18)+" $Å^3$")
        if unit_cell_volume != 0:
            with col1:
                st.write("Unit cell ÷ Hofmann =", str(round(
                                        unit_cell_volume / total_hofmann, 3)))
            #with col2:
            #    st.write("Unit cell ÷ 18Å =", str(round(unit_cell_volume / total_18, 2)))
elif option == "Examples":
    with st.beta_expander(label="Example inputs", expanded=True):
        """
        There are a number of acceptable inputs that can be used to estimate
        crystallographic volumes using Hofcalc.

        The simplest option is to enter the chemical formula, name or SMILES
        of the material of interest, e.g.
        """
        search_terms = [["CH3CH2OH", "Formula", "69.61"],
                        ["ethanol", "name", "69.61"], ["OCC", "SMILES", "69.61"]]
        df = pd.DataFrame(search_terms, columns=["Search term","Type","Volume"])
        st.table(df)
        st.write("")
        """
        It is also possible to search for multiple items of any type at
        the same time, provided that they can be separated by a space e.g.
        """
        search_terms = [["carbamazepine indomethacin", "731.97"],
                        ["zopiclone 2H2O", "496.02"], ["C15H12N2O CH3CH2COO- Na+",
                        "419.79"]]
        df = pd.DataFrame(search_terms, columns=["Search term", "Volume"])
        st.table(df)
        """
        However, this can cause problems when the name of a compound requires
        multiple words, for example "acetic acid".

        If this occurs, uncheck the  "Autodetect Fragments" checkbox, and
        try again.

        It is possible to manually add additional fragments using the
        Number of fragments element that appears when the checkbox is unchecked
        """
        search_terms = [["acetic acid", "70.84"],
                        ["verapamil hydrochloride", "667.57"], ["sodium salicylate",
                        "182.66"]]
        df = pd.DataFrame(search_terms, columns=["Search term", "Volume"])
        st.table(df)
    for i in range(3):
        st.write("")
    with st.beta_expander(label="Temperature", expanded=True):
        """
        The temperature (in kelvin) can be entered, which will automatically be
        applied to the volume calculation using the follwing equation:
        """
        st.latex("V = \\sum{n_{i}v_{i}(1 +  \\alpha(T - 298))}")
        """
        Where
        """
        st.latex("\\alpha = 0.95 \\times 10^{-4}")
    for i in range(3):
        st.write("")
    with st.beta_expander(label="Unit Cell", expanded=True):
        """
        If the volume of a unit cell is supplied, then in addition to displaying
        the estimated volume, the cell volume divided by the estimated molecular
        volume will also be shown.
        """
        search_terms = [["zopiclone 2H2O", "1874.61", "496.02", "3.78"],
                        ["verapamil HCl", "1382.06", "667.57", "2.07"]]
        df = pd.DataFrame(search_terms, columns=["Search term",
                    "Unit Cell Volume", "Hofmann Volume", "Vcell / VHofmann"])
        st.table(df)
else:
    st.write("Average crystallographic volumes at 298 K reported by Hofmann \
            (see references)")
    volume_df = pd.DataFrame.from_dict(volumes, columns=["Volume"],
                                        orient="index", dtype=float)
    st.table(volume_df.style.format('{:.2f}'))


