import streamlit as st
from pyvalem.formula import Formula
from collections import defaultdict
import pubchempy as pcp
import math
import json

with open("volumes.json", "r") as f:
    volumes = json.load(f)
f.close()

def get_density(formula, temperature=298):
    mass = 0
    for element in formula.keys():
        mass += Formula(element).rmm * formula[element]
    mass *= 1.66054e-24 # convert mass to grams
    volume = get_volume(formula, temperature=temperature)
    volume *= (1e-8)**3 # convert volume to cm^3
    return round(mass / volume, 2)

def get_volume(formula, temperature=298):
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
                    # Use this to trip up deuterium, tritium issues when
                    # searching for common abbreviations e.g. THC, CBD, which
                    # could be interpreted as chemical formulae, but don't have
                    # corresponding Hofmann tabulated values (identical to H)
                    for element in mf.keys():
                        _ = volumes[element]
                except:
                    try:
                        component_no_decimal = component.replace(".","1",1)
                        i = 0
                        for c in component_no_decimal:
                            if c.isdigit():
                                i+=1
                            else:
                                break
                        multiple = float(component[:i])
                        actual_formula = component[i:]
                        f = Formula(actual_formula)
                        mf = f.atom_stoich
                        # Use this to trip up deuterium, tritium issues when
                        # searching for common abbreviations e.g. THC, CBD, which
                        # could be interpreted as chemical formulae, but don't have
                        # corresponding Hofmann tabulated values (identical to H)
                        for element in mf.keys():
                            _ = volumes[element]
                        for key, value in mf.items():
                            mf[key] = value * multiple

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