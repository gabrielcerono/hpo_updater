import requests


def request_get(url, headers):
    """
    Performs a get request that provides a (somewhat) useful error message.
    """
    try:
        response = requests.get(url, headers=headers)
    except ImportError:
        raise ImportError("Couldn't retrieve the data, check your URL")
    else:
        return response


def get_json(url, headers):
    """Returns request in JSON (dict) format"""
    response = request_get(url, headers)
    if response.status_code != 204 and\
       response.headers["Content-Type"].strip().startswith("application/json"):
        return response.json()
    else:
        raise Exception("Response is not a JSON object")
