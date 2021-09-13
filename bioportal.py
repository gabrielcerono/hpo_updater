import os
from bioportal_client import BioPortalClient


def get_bioportal_mappings(from_ontology, to_ontology):
    api_key = os.getenv('BIOPORTAL_API_KEY')
    client = BioPortalClient(api_key)
    return client.get_mappings(from_ontology, to_ontology)
