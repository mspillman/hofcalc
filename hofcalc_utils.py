import streamlit as st
from pyvalem.formula import Formula
from collections import defaultdict
import pubchempy as pcp
import math

volumes = {
    "H" : 5.08, "He" : math.nan, "Li" : 22.6, "Be" : 36., "B" : 13.24,
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
    "Es" : math.nan, "Fm" : math.nan
    }

def get_density(formula, temperature=298):
    mass = 0
    for element in formula.keys():
        mass += Formula(element).rmm * formula[element]
    mass *= 1.66054e-24
    volume = get_volume(formula, temperature=temperature)
    volume *= (1e-8)**3
    return round(mass / volume, 2)


def get_volume(formula, volumes=volumes, temperature=298):
    volume = 0.0
    alpha = 0.95e-4
    eighteen_angstrom_volume = 0
    for element in formula.keys():
        if element not in volumes:
            st.write("Element",element,"not recognized.")
            return None
        else:
            if element != "H":
                eighteen_angstrom_volume += 18*formula[element]
            if volumes[element] is not math.nan:
                volume += (volumes[element] *
                            formula[element])*(1.+alpha*(temperature - 298.))
            else:
                st.write("Element", element, "does not have a Hofmann volume")
                st.write("Using the 18A rule to estimate missing volume")
                volume += 18. * formula[element]
    return round(volume, 2)

def get_formula():
    user_input = st.text_input("Enter chemical formula or name. Multiple \
                    entries should be separated by semicolons.",
                    value='', max_chars=None, key=None, type='default')
    if user_input != "":
        cids = {}
        cid_links = {}
        # Remove html tags produced by MarvinSketch
        if "<" in user_input:
            user_input = user_input.split("</font>")[0]
            user_input = user_input.split("<b>")[1]
            user_input = user_input.split(">")[1:]
            user_input = ("".join(user_input)).replace("<sub","").replace(
                                                    "</sub","").replace(" ","")
            st.write("Interpreting formula as", user_input)
        molecular_formulae = {}
        if ";" in user_input:
            components = user_input.split(";")
            components = [x.strip() for x in components]
        else:
            components = [user_input.strip()]
        for component in components:
            if component != "":
                try:
                    component = component.replace("|", "")
                    f = Formula(component)
                    mf = f.atom_stoich
                except:
                    try:
                        molecule = pcp.get_compounds(component, "name")[0]
                        cids[component] = molecule.cid
                        base_url = "https://pubchem.ncbi.nlm.nih.gov/compound/"
                        cid_links[component] = base_url+str(molecule.cid)
                        f = Formula(molecule.molecular_formula)
                        mf = f.atom_stoich
                    except:
                        st.write("Error")
                        st.write("Unable to parse",component,"as a chemical \
                            formula or find it on PubChem when searching by \
                            name.")
                        return None
                molecular_formulae[component] = mf

        combined_molecular_formula = defaultdict(int)
        for component in components:
            for key, value in molecular_formulae[component].items():
                combined_molecular_formula[key] += value
        molecular_formula = {"combined" : combined_molecular_formula,
                            "individual" : molecular_formulae,
                            "user_input" : components,
                            "PubChem CIDs" : cids,
                            "PubChem URLs" : cid_links}
    else:
        molecular_formula = None

    return molecular_formula