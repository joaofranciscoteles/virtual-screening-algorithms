from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import scaffoldgraph as sg
import pandas as pd

def extract_scaffold(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        return None

    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)

    if scaffold.GetNumAtoms() == 0:
        return None

    return Chem.MolToSmiles(scaffold)

def extract_scaffolds(molecules: list[dict]) -> list[dict]:
    results = []

    for item in molecules:
        scaffold = extract_scaffold(item["smiles"])

        results.append({
            "id": item["id"],
            "smiles": item["smiles"],
            "scaffold": scaffold,
        })

    return results


def build_hierarchy(molecules: list[dict]):
    dataframe = pd.DataFrame(molecules)

    graph = sg.ScaffoldNetwork.from_dataframe(
        dataframe,
        smiles_column="smiles",
        name_column="id",
    )

    return graph

def analyze_scaffolds(molecules: list[dict]):
    extracted = extract_scaffolds(molecules)
    hierarchy = build_hierarchy(molecules)

    return {
        "molecules": extracted,
        "hierarchy": hierarchy,
    }