from rdkit import Chem

def substructure_search(molecule_smiles_list, substructure_smiles):
    """
   Search for molecules in smiles_list that contain the substructure specified by substructure_smiles.

   Args:
   smiles_list (list of str): List of SMILES strings representing molecules.
   substructure_smiles (str): SMILES string representing the substructure.

   Returns:
   list of str: List of SMILES strings of molecules that contain the substructure.
    """
    substructure = Chem.MolFromSmiles(substructure_smiles)
    matching_molecules = []

    for smiles in molecule_smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule.HasSubstructMatch(substructure):
            matching_molecules.append(smiles)

    return matching_molecules


if __name__ == "__main__":
    molecule_smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure_smiles = "c1ccccc1"
    result = substructure_search(molecule_smiles_list, substructure_smiles)
    print(result)
    print("hii")
