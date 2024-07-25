from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

app = FastAPI()


molecules_db = {}



class Molecule(BaseModel):
    identifier: str
    smiles: str


class SubstructureSearch(BaseModel):
    query_smiles: str


@app.get("/molecules/{identifier}")
def get_molecule(identifier: str):
    smiles = molecules_db.get(identifier)
    if smiles is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": identifier, "smiles": smiles}

@app.post("/molecules")

def add_molecule(molecule: Molecule):
    if molecule.identifier in molecules_db:
        raise HTTPException(status_code=400, detail="Molecule already exists")
    molecules_db[molecule.identifier] = molecule.smiles
    return {"message": "Molecule added successfully"}


@app.put("/molecules/{identifier}")
def update_molecule(identifier: str, molecule: Molecule):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecules_db[identifier] = molecule.smiles
    return {"message": "Molecule updated successfully"}


@app.delete("/molecules/{identifier}")
def delete_molecule(identifier: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules_db[identifier]
    return {"message": "Molecule deleted successfully"}


@app.get("/molecules")
def list_molecules():
    return molecules_db


@app.post("/substructure_search")
def substructure_search(search: SubstructureSearch):
    query = Chem.MolFromSmiles(search.query_smiles)
    if query is None:
        raise HTTPException(status_code=400, detail="Invalid query SMILES")

    matches = {}
    for identifier, smiles in molecules_db.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(query):
            matches[identifier] = smiles
    return matches



@app.post("/upload")
def upload_file(file: bytes):
    # Implement file parsing and molecule addition here
    return {"message": "File uploaded successfully"}
