import streamlit as st
import pandas as pd
import hofcalc_utils as hu
import json
import base64

st.set_page_config(page_title='HofCalc WebApp', page_icon = ":test_tube:",
                    layout = 'centered', initial_sidebar_state = 'auto')

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.title("HofCalc")
st.markdown("*Hofmann Volume Calculator*")
st.sidebar.title("Options")

with st.sidebar:
    option = st.radio("", ["Volume Estimation",
                            "Hofmann Volumes (298 K)",
                            "Help"])

if option == "Volume Estimation":
    col1, col2 = st.beta_columns(2)

    with col1:
        temperature = st.number_input("Temperature / K",
                                min_value=0, max_value=None, value=298, step=10)
    with col2:
        unit_cell_volume = st.number_input(
                                "Unit cell volume in cubic ångströms (optional)",
                                min_value=0.0, max_value=None, value=0.0,
                                step=100.0)
    molecular_formula = hu.get_formula()
    total_volume = 0
    st.write("")
    st.write("")
    if molecular_formula is not None:
        st.markdown("**Results**")
        if len(molecular_formula["individual"].keys()) > 1:
            with st.beta_expander(label="Individual", expanded=False):
                individual_volumes = {}
                col1, col2, col3 = st.beta_columns(3)
                with col1:
                    st.markdown("***Input***")
                with col2:
                    st.markdown("***Formula***")
                with col3:
                    st.markdown("***Volume***")
                for mf in molecular_formula["individual"].keys():
                    individual_volumes[mf] = hu.get_volume(
                                            molecular_formula["individual"][mf],
                                            temperature=temperature)
                for mf, volume in zip(molecular_formula["individual"],
                                    individual_volumes.values()):
                    with col1:
                        st.write(mf)
                    with col2:
                        mf = molecular_formula["individual"][mf]
                        mf_string = ""
                        for key, val in mf.items():
                            mf_string += key + str(val) + " "
                        st.write(mf_string)
                    with col3:
                        st.write(str(volume))
            for i in range(3):
                st.write("")
            molecular_formula["individual_volumes"] = individual_volumes
        total_volume = hu.get_volume(molecular_formula["combined"],
                                        temperature=temperature)

        col1, col2, col3, col4 = st.beta_columns(4)
        with col1:
            mf = molecular_formula["combined"]
            mf_string = ""
            for key, val in mf.items():
                mf_string += key + str(val) + " "
            st.write("Total atoms:")
            st.write(mf_string)
        with col2:
            st.write("Hofmann Volume:")
            st.markdown(str(round(total_volume, 3))+" $Å^3$")
        with col3:
            density = hu.get_density(molecular_formula["combined"], temperature)
            st.write("Hofmann Density:")
            st.write(str(density),"$g$ $cm^{-3}$")
        if unit_cell_volume != 0:
            with col4:
                st.write("$$\\frac{V_{Cell}}{V_{Hofmann}}$$:")
                st.write(str(round(unit_cell_volume / total_volume, 2)))
                molecular_formula["V_Cell / V_Hofmann"] = round(
                                        unit_cell_volume / total_volume, 2)
        for i in range(2):
            st.write("")

        molecular_formula["Temperature"] = temperature
        molecular_formula["Hofmann Volume"] = total_volume
        molecular_formula["Hofmann Density"] = density

        name = " ".join([x.strip() for x in molecular_formula["user_input"]])
        name = name.replace(" ","_")
        name += "_"+str(round(temperature, 2))+"K_HofCalc.json"
        formula_json = json.dumps(molecular_formula, indent=4)
        b64 = base64.b64encode(formula_json.encode()).decode()
        href = f'<a href="data:file/zip;base64,{b64}" download=\'{name}\'>\
                Download summary</a>'
        st.markdown(href, unsafe_allow_html=True)

elif option == "Help":
    st.markdown("A blog post describing HofCalc and its use can be found [here](https://mspillman.github.io/blog/pxrd/indexing/2021/11/10/Hofcalc.html)")
    with st.beta_expander(label="Example inputs", expanded=False):
        st.markdown("**Basic**")
        """
        There are a number of acceptable inputs that can be used to estimate
        crystallographic volumes using HofCalc.

        The simplest option is to enter the chemical formula or name of the
        material of interest. Names are resolved by querying PubChem.
        Note that formulae can be prefixed with a multiple, e.g. 2H2O
        """
        search_terms = [["ethanol", "name", "69.61"],
                        ["CH3CH2OH", "formula", "69.61"],
                        ["water", "name", "21.55"],
                        ["2H2O", "formula", "43.10"]]
        df = pd.DataFrame(search_terms, columns=["Search term","Type","Volume"])
        st.table(df)
        st.write("")
        st.markdown("**Multiple search terms**")
        """
        It is also possible to search for multiple items simultaneously, and mix
        and match name and formulae by separating individual components with a
        semicolon. This means that for example, 'amodiaquine dihydrochloride
        dihydrate' can also be entered as 'amodiaquine; 2HCl; 2H2O'.
        """
        search_terms = [["carbamazepine; L-glutamic acid", "497.98"],
                        ["zopiclone; 2H2O", "496.02"],
                        ["C15H12N2O; CH3CH2COO-; Na+", "419.79"],
                        ["sodium salicylate; water", "204.21"],
                        ["amodiaquine dihydrochloride dihydrate", "566.61"],
                        ["amodiaquine; 2HCl; 2H2O", "566.61"]]
        df = pd.DataFrame(search_terms, columns=["Search term", "Total Volume"])
        st.table(df)
        st.markdown("**More complex examples - hemihydrates**")
        """
        In cases where fractional multiples of search components are required,
        such as with hemihydrates, care should be taken to check the evaluated
        chemical formula for consistency with the expected formula.
        """
        search_terms = [["Calcium sulfate hemihydrate","Ca2 H2 O9 S2","253.07",
                        "2", "126.53"],
                        ["calcium; calcium; sulfate; sulfate; water",
                        "Ca2 H2 O9 S2", "253.07", "2", "126.53"],
                        ["calcium; sulfate; 0.5H2O", "Ca1 H1.0 O4.5 S1",
                        "126.53", "-", "126.53"],
                        ["Codeine phosphate hemihydrate", "C36 H50 N2 O15 P2",
                        "1006.77", "2", "503.38"],
                        ["codeine; codeine; phosphoric acid; \
                        phosphoric acid; water", "C36 H50 N2 O15 P2",
                        "1006.77", "2", "503.38"],
                        ["codeine; phosphoric acid; 0.5H2O",
                        "C18 H25.0 N1 O7.5 P1", "503.38", "-", "503.38"]]
        df = pd.DataFrame(search_terms, columns=["Search term", "Evaluated as",
                        "Calculated Volume", "Divide by", "Expected Volume"])
        st.table(df)
        st.markdown("**Charged species in formulae**")
        """
        Charges could potentially interfere with the parsing of chemical
        formulae. For example, two ways of representing an oxide ion:
        """
        search_terms = [["O-2", "1 x O"],
                        ["O2-", "2 x O"]]
        df = pd.DataFrame(search_terms, columns=["Search term", "Evaluated as"])
        st.table(df)
        """
        Whilst is is recommended that charges be omitted from HofCalc queries,
        if including charges in your queries, ensure that the correct number of
        atoms has been determined in the displayed atom counts or the
        downloadable summary file. For more information on formatting formulae,
        see the pyvalem documentation (link in references)
        """
    for i in range(3):
        st.write("")
    with st.beta_expander(label="Temperature", expanded=False):
        """
        The temperature, $T$ (in kelvin) is automatically included in the volume
        calculation via the following equation:
        """
        st.latex("V = \\sum{n_{i}v_{i}(1 +  \\alpha(T - 298))}")
        """
        Where $n_{i}$ and $v_{i}$ are the number and Hofmann volume (at 298 K)
        of the $i$th element in the chemical formula, and where
        """
        st.latex("\\alpha = 0.95 \\times 10^{-4} K^{-1}")

    for i in range(3):
        st.write("")
    with st.beta_expander(label="Unit Cell", expanded=False):
        """
        If the volume of a unit cell is supplied, then the unit cell volume
        divided by the estimated molecular volume will also be shown.
        """
        search_terms = [["zopiclone, 2H2O", "1874.61", "496.02", "3.78"],
                        ["verapamil, HCl", "1382.06", "667.57", "2.07"]]
        df = pd.DataFrame(search_terms, columns=["Search term",
                    "Unit Cell Volume", "Hofmann Volume", "Vcell / VHofmann"])
        st.table(df)
    for i in range(3):
        st.write("")
    with st.beta_expander(label="Summary files", expanded=False):
        """
        Each time HofCalc is used, a downloadable summary file is produced.
        It is designed to serve both as a record of the query for future
        reference and also as a method to sense-check the interpretation of the
        entered terms, with links to the PubChem entries where relevant.

        An example of the contents of the summary file for the following search
        terms is given below.

        "carbamazepine; indomethacin" with T = 293 and unit cell volume = 2921.6
        """
        summary = {
            "combined": {
                "C": 34,
                "H": 28,
                "N": 3,
                "O": 5,
                "Cl": 1
            },
            "individual": {
                "carbamazepine": {
                    "C": 15,
                    "H": 12,
                    "N": 2,
                    "O": 1
                },
                "indomethacin": {
                    "C": 19,
                    "H": 16,
                    "Cl": 1,
                    "N": 1,
                    "O": 4
                }
            },
            "user_input": [
                "carbamazepine",
                "indomethacin"
            ],
            "PubChem CIDs": {
                "carbamazepine": 2554,
                "indomethacin": 3715
            },
            "PubChem URLs": {
                "carbamazepine": "https://pubchem.ncbi.nlm.nih.gov/compound/2554",
                "indomethacin": "https://pubchem.ncbi.nlm.nih.gov/compound/3715"
            },
            "individual_volumes": {
                "carbamazepine": 303.86,
                "indomethacin": 427.77
            },
            "V_Cell / V_Hofmann": 3.99,
            "Temperature": 293,
            "Hofmann Volume": 731.62,
            "Hofmann Density": 1.35
        }
        st.write(summary)
    for i in range(3):
        st.write("")
    with st.beta_expander(label="References", expanded=False):
        col1, col2 = st.beta_columns([2,2])
        with col1:
            st.write("Hofmann, D.W.M. (2002), Fast estimation of crystal \
                    densities. Acta Cryst. B, 58: 489-493.")
        with col2:
            st.write("https://doi.org/10.1107/S0108768101021814")
        st.write("")
        st.write("")
        col1, col2 = st.beta_columns([2,2])
        with col1:
            st.write("PubChem")
        with col2:
            st.write("https://pubchem.ncbi.nlm.nih.gov/")
        st.write("")
        st.write("")
        col1, col2 = st.beta_columns([2,2])
        with col1:
            st.write("PubChemPy")
        with col2:
            st.write("https://github.com/mcs07/PubChemPy")
        st.write("")
        st.write("")
        col1, col2 = st.beta_columns([2,2])
        with col1:
            st.write("PyValem")
        with col2:
            st.write("https://github.com/xnx/pyvalem")
        st.write("")
        st.write("")
        col1, col2 = st.beta_columns([2,2])
        with col1:
            st.write("HofCalc WebApp")
        with col2:
            st.write("https://github.com/mspillman/hofcalc")

else:
    st.write("Average crystallographic volumes at 298 K reported by Hofmann \
            (see references)")
    st.write("All volumes are in cubic ångströms")
    volume_df = pd.DataFrame.from_dict(hu.volumes, columns=["Volume"],
                                        orient="index", dtype=str)
    st.table(volume_df)


