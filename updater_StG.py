import pandas as pd
import requests
import os
import subprocess
from bioportal_client import BioPortalClient


def get_bioportal_mappings(from_ontology, to_ontology):
    api_key = "225a7f07-744e-41fc-b8e8-7c3d9b11a100"
    client = BioPortalClient(api_key)
    return client.get_mappings(from_ontology, to_ontology)

def robot_input():
    with open("get_dbxref.txt", 'w') as f:
	    f.write('PREFIX owl: <http://www.w3.org/2002/07/owl#>\nPREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\nPREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\nPREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>\n\nSELECT ?concept ?label ?dbxref WHERE {\n	?concept rdf:type owl:Class ;\n	         rdfs:label ?label ;\n		     oboInOwl:hasDbXref ?dbxref .\n    FILTER(regex(?dbxref, "MSH:"))\n}')

def robot_bash():
    """ Runs the robot program to return the MPO-MESH mappings"""
    robot_input()
    subprocess.run('./robot convert -I http://purl.obolibrary.org/obo/hp.owl -o hp.owl extract -b http://purl.obolibrary.org/obo/HP_0000118 -m MIREOT query --query get_dbxref.txt get_dbxref.csv', shell=True)

def hpo_mesh_map():
    """ Get Hpo to Mesh mappings and parse into a Dataframe from ROBOT """

    robot_bash()
    
    hpo_mesh_df = pd.read_csv("get_dbxref.csv")
    hpo_mesh_df['HPO']=hpo_mesh_df['concept'].apply(lambda row: row.split('_')[-1])
    hpo_mesh_df['MESH']=hpo_mesh_df['dbxref'].apply(lambda row: row.split(':')[-1])
    hpo_mesh_df.drop(columns = ['concept', 'dbxref'], inplace= True)
    hpo_mesh_df.rename(columns={'label' : "mesh_name"}, inplace= True)
    return hpo_mesh_df

def get_hpo_annot():
    """ Get HPO annotations from HPO jax org, these are edges from
    HPO terms to Orpha and Omim """
    url = 'http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'
    hpo_annot = pd.read_csv(url, delimiter = "\t", skiprows = 2, header = None)
    header_key = ["HPO_ID", "HPO_LABEL", "gen_id", "gen_name", "info_from_source", "source", "disease_id"]
    hpo_annot.columns = header_key
    hpo_annot.drop(columns = 'info_from_source', inplace = True)
    return hpo_annot

def get_edges():
    """ Get edges from Mesh to genes. We get this by merging on HPO annotations and HPO/MESH mappings. """

    print('getting hpo annotations')
    hpo_annot = get_hpo_annot()
    print('getting hpo - mesh mappings')
    hpo_mesh_df = hpo_mesh_map()
    print('generating the mesh-genes edges')
    hpo_annot['HPO'] = hpo_annot['HPO_ID'].apply(lambda row : row.split(':')[-1])
    hpo_annot.drop(columns= 'HPO_ID', inplace= True)
    hpo_mesh = hpo_annot.merge(hpo_mesh_df, how= "inner", on = 'HPO')
    hpo_mesh.drop(columns = ['HPO', 'disease_id', 'HPO_LABEL', 'gen_name'], inplace= True)
    hpo_mesh.drop_duplicates(inplace = True)
    return hpo_mesh

def delete_bash():
    subprocess.run("rm -f get_dbxref.txt get_dbxref.csv hp.owl", shell= True)


edges = get_edges()
edges.to_csv('edges.csv')
delete_bash()
