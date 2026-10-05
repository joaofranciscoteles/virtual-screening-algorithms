#fingerprint --> transformar a molécula em números
#maccs e morgan/ecfp sao duas maneiras de fazer isso
#tanimoto --> aqui é para comparar dois fingerprints 
# no final vou rankear os resultados 

#Dúvidas para confirmar: no artigo tabela 2, a entrada está como binary/count vector, mas isso seria a entrada para o tanimoto,certo?pois o morgan e maccs sao os reponsaveis por gerar os fingerprints, que são vetores binários ou de contagem, e o Tanimoto é a métrica que compara esses vetores para calcular a similaridade entre moléculas. 
#aqui eu usei como entrada esse molecule, que é uma estrutura ja interpretada do rdkit, e essa estrutura pode vir de smiles ou de outros formatos, preciso confirmar se é isso mesmo, se sera essa entrada mesmo? 

from rdkit import DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

def generate_morgan(molecule, radius, size):

    if molecule is None:
        raise ValueError("Invalid molecule.")

    if radius < 0 or size <= 0:
        raise ValueError("Radius must be >= 0 and size must be > 0.")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=size,
    )

    return generator.GetFingerprint(molecule)


def generate_maccs(molecule):
    
    if molecule is None:
        raise ValueError("Invalid molecule.")

    return MACCSkeys.GenMACCSKeys(molecule)


def calculate_tanimoto(reference_fp, candidate_fp):
    
    return DataStructs.TanimotoSimilarity(
        reference_fp,
        candidate_fp,
    )


def generate_ranking(reference_fp, candidates):
    
    results = []

    for identifier, candidate_fp in candidates:
        similarity = calculate_tanimoto(
            reference_fp,
            candidate_fp,
        )

        results.append({
            "id": identifier,
            "similarity": similarity,
        })

    results.sort(
        key=lambda result: result["similarity"],
        reverse=True,
    )

    return results