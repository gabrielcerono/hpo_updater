import os 
from bioportal_client import BioPortalClient
import pandas as pd
import requests


def get_bioportal_mappings(from_ontology, to_ontology):
    api_key = "225a7f07-744e-41fc-b8e8-7c3d9b11a100"
    client = BioPortalClient(api_key)
    return client.get_mappings(from_ontology, to_ontology)


def hpo_mesh_map():
    """ Get Hpo to Mesh mappings and parse into a Dataframe """

    hpo_mesh = get_bioportal_mappings('HP', 'MESH')
    hpo_mesh_df = pd.DataFrame.from_dict(hpo_mesh, orient='index', columns = ['HPO'])
    hpo_mesh_df.reset_index(inplace = True)
    hpo_mesh_df['HPO']=hpo_mesh_df['HPO'].apply(lambda row: row.split('_')[-1])
    hpo_mesh_df['MESH']=hpo_mesh_df['index'].apply(lambda row: row.split('/')[-1])
    hpo_mesh_df.drop(columns = 'index', inplace = True)
    return hpo_mesh_df

def omim_do_map():
    """ Get Omim to Doid mappings and parse into a Dataframe """

    doid_omim = get_bioportal_mappings("DOID", "OMIM")
    do_omim_df = pd.DataFrame.from_dict(doid_omim, orient='index', columns = ['OMIM'])
    do_omim_df.reset_index(inplace = True)
    do_omim_df['DOID']=do_omim_df['index'].apply(lambda row: row.split('_')[-1])
    do_omim_df['OMIM']=do_omim_df['OMIM'].apply(lambda row: row.split('/')[-1])
    do_omim_df.drop(columns = 'index', inplace = True)
    return do_omim_df

def orpha_do_map():
    """ Get Orpha to Doid mappings and parse into a Dataframe"""
    doid_orpha = get_bioportal_mappings('DOID', 'ORDO')
    do_orpha_df = pd.DataFrame.from_dict(doid_orpha, orient='index', columns = ['ORPHA'])
    do_orpha_df.reset_index(inplace = True)
    do_orpha_df['DOID']=do_orpha_df['index'].apply(lambda row: row.split('_')[-1])
    do_orpha_df['ORPHA']=do_orpha_df['ORPHA'].apply(lambda row: row.split('_')[-1])
    do_orpha_df.drop(columns = 'index', inplace = True)
    return do_orpha_df

def get_hpo_annot():
    """ Get HPO annotations from HPO jax org, these are edges from
    HPO terms to Orpha and Omim """
    url = 'http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'
    hpo_annot = pd.read_csv(url, delimiter = "\t", skiprows = 2, header = None)
    header_key = ["HPO_ID", "HPO_LABEL", "gen_id", "gen_name", "info_from_source", "source", "disease_id"]
    hpo_annot.columns = header_key
    hpo_annot.drop(columns = 'info_from_source', inplace = True)
    return hpo_annot



def main():
    """ Build a Dataframe with 2 Columns, one is the Mesh term and the Other is the DOID"""
    #Get the edges 
    hpo_annot = get_hpo_annot()

    #Get only the edges to OMIM
    hpo_edge_omim = hpo_annot.loc[hpo_annot['disease_id'].str.contains('OMIM', na = False)]
    hpo_edge_omim['OMIM'] = hpo_edge_omim['disease_id'].apply(lambda row : row.split(':')[-1])
    hpo_edge_omim['HPO'] = hpo_edge_omim['HPO_ID'].apply(lambda row : row.split(':')[-1])

    #Merge edges from HPO to OMIM, with HPO-MESH mappings
    hpo_mesh_df = hpo_mesh_map()
    hpo_mesh_edge_to_omim = hpo_edge_omim.merge(hpo_mesh_df, how = 'inner', on = 'HPO')
    hpo_mesh_edge_to_omim['OMIM'] = hpo_mesh_edge_to_omim['disease_id'].apply(lambda row : row.split(':')[-1])

    #We need another merge, with OMIM-DO mappings (Right now we have hpo/mesh to OMIM)
    do_omim_df = omim_do_map()
    hpo_mesh_edge_to_omim_do = hpo_mesh_edge_to_omim.merge(do_omim_df, how = 'inner', on = 'OMIM')

    #Clean the Dataframe dropping duplicates and unnecesary columns.
    hpo_mesh_edge_to_omim_do.drop(columns=['gen_name', 'gen_id', 'HPO_ID', 'source', 'disease_id', "HPO_LABEL", "OMIM", "HPO"], inplace = True)
    hpo_mesh_edge_to_omim_do.drop_duplicates(inplace=True)

    #Orpha part of the HPO annotations
    hpo_edge_orpha = hpo_annot.loc[hpo_annot['disease_id'].str.contains('ORPHA', na = False)]
    hpo_edge_orpha['ORPHA'] = hpo_edge_orpha['disease_id'].apply(lambda row : row.split(':')[-1])
    hpo_edge_orpha['HPO'] = hpo_edge_orpha['HPO_ID'].apply(lambda row : row.split(':')[-1])

    #Merge edges from HPO to Orpha, with HPO-MESH mappings
    hpo_mesh_edge_to_orpha = hpo_edge_orpha.merge(hpo_mesh_df, how = 'inner', on = 'HPO')

    #Merge edges from HPO/MESH to Orpha, with ORPHA/DOID mappings
    do_orpha_df = orpha_do_map()
    hpo_mesh_edge_to_orpha_do = hpo_mesh_edge_to_orpha.merge(do_orpha_df, how = 'inner', on = 'ORPHA')

    #Clean the Dataframe dropping unnecesary columns and duplicates
    hpo_mesh_edge_to_orpha_do.drop(columns = ['gen_name', 'gen_id', 'HPO_ID', 'source', 'disease_id', 'HPO_LABEL', 'ORPHA', 'HPO'], inplace = True)
    hpo_mesh_edge_to_orpha_do.drop_duplicates(inplace = True)

    mesh_edge_to_do = pd.concat([hpo_mesh_edge_to_orpha_do, hpo_mesh_edge_to_omim_do], ignore_index= True)

    return mesh_edge_to_do



