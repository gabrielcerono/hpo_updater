from request_utils import get_json


class BioPortalClient:
    """BioPortal API client
    Provides functions to easily access the BioPortal API
    (http://data.bioontology.org/documentation) in Python.

    Attributes:
        get_mappings: Retrieves BioPortal mappings given two ontologies.
    """

    _BASE_URL = "https://data.bioontology.org"
    _MAPPINGS_API = "mappings"

    _ONTOLOGIES = "ontologies"
    _PAGE_SIZE = "pagesize"

    def __init__(self, api_key):
        self.api_key = api_key
        self.headers = {'Authorization': 'apiKey token=' + api_key}

    def get_mappings(self, from_ontology, to_ontology):
        """Returns a Python dictionary that contains a mapping between
           concepts in the from_ontology and to_ontology.

        Args:
            from_ontology (str): The source ontology identifier from BioPortal.
            to_ontology (str): The target ontology identifier from BioPortal

        Returns:
            A Python dictionary.
        """
        mappings_api_base = f"{self._BASE_URL}/{self._MAPPINGS_API}"
        ontologies_option = f"{self._ONTOLOGIES}={to_ontology},{from_ontology}"
        page_size_option = f"{self._PAGE_SIZE}=100"
        url = f"{mappings_api_base}?{ontologies_option}&{page_size_option}"

        mappings = {}

        json = get_json(url, self.headers)
        while json['links']['nextPage'] is not None:
            mappings.update(self._parse_mappings(json))
            next_page = json['links']['nextPage']
            json = get_json(next_page, self.headers)

        mappings.update(self._parse_mappings(json))

        return mappings

    def _parse_mappings(self, json):
        mapping = {}
        collection = json['collection']
        for item in collection:
            map_item = item['classes']
            key = map_item[0]['@id']
            value = map_item[1]['@id']
            mapping[key] = value
        return mapping
