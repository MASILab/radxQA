#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Utility to sync kurts stuff to a redcap database """

# pylint: disable-msg=too-many-locals, too-many-statements

import os
from argparse import ArgumentParser

import requests

AP = ArgumentParser(description='Utility to sync kurts stuff to a redcap database',
                    epilog='Contact: kurt.g.schilling.1@vanderbilt.edu')
AP.add_argument('csv_path')
AP.add_argument('series_image')

def get_records(api_url, api_key):
    """ returns records of a redcap project given an api url and api key """
    payload = {'token': api_key,
               'content': 'record',
               'format': 'csv',
               'type': 'flat'}

    response = requests.post(api_url, data=payload)

    if response.status_code == 200:
        return response.content.decode(response.encoding)
    else:
        # status codes are incredibly confusing and don't help debugging in my experience
        # so just remind user to check api url, api key, user permissions and host
        raise RuntimeError('Could not get records. Please verify that the API URL and API KEY \
                            are correct, proper user permissions are set, and that the host is up.')

def set_records(api_url, api_key, data):
    """ set records of a redcap project given an api url, api key, and data """
    payload = {'token': api_key,
               'content': 'record',
               'format': 'csv',
               'type': 'flat',
               'overwriteBehavior': 'overwrite',
               'forceAutoNumber': 'true',
               'returnContent': 'auto_ids',
               'data': data}

    # If forceAutoNumber is set to true, redcap will automatically assign a record number.
    # However, record numbers during the import are still required to associate rows to the same
    # records

    response = requests.post(api_url, data=payload)

    if response.status_code == 200:
        return response.content.decode(response.encoding)
    else:
        # status codes are incredibly confusing and don't help debugging in my experience
        # so just remind user to check api url, api key, user permissions and host
        raise RuntimeError('Could not set records. Please verify that the API URL and \
                            API KEY are correct, proper user permissions are set, and that the \
                            host is up.')

def set_file(api_url, api_key, file_path, record, field):
    """ sets a file of a redcap project given an api url, api key, file path, record, and field """
    with open(file_path, 'rb') as file:
        payload = {'token': api_key,
                   'content': 'file',
                   'action': 'import',
                   'record': record,
                   'field': field}
        files = {'file': (file_path, file)}

        response = requests.post(api_url, data=payload, files=files)

        if response.status_code == 200:
            return response.content.decode(response.encoding)
        else:
            # status codes are incredibly confusing and don't help debugging in my experience
            # so just remind user to check api url, api key, user permissions and host
            raise RuntimeError('Could not set file. Please verify that the API URL and \
                                API KEY are correct, proper user permissions are set, and that the \
                                host is up.')

def main():
    """ Main program """

    # Get API url and API key
    #api_url = os.environ.get('API_URL')
    #if not api_url:
    #    raise RuntimeError('"API_URL" was not set as an environment variable.')
    api_url = 'https://redcap.vanderbilt.edu/api/'

    #api_key = os.environ.get('API_KEY')
    #if not api_key:
    #    raise RuntimeError('"API_KEY" was not set as an environment variable.')
    api_key = '7597F8320EB167576D9C39381768A633'
    
    # Handle arguments
    args = AP.parse_args()

    csv_path = args.csv_path
    series_image = args.series_image

    # Upload record
    print('Uploading record...')
    with open(csv_path, 'rb') as file:
        response = set_records(api_url, api_key, file.read())

    # Get record
    #records = get_records(api_url, api_key)
    #record = records.split('\n')[-2].split(',')[0]
    response_ids = [response_split.split(',')[0] for response_split in response.split('\n')[1:]]
    print('Record #: ' + str(response_ids))

    # Upload series_image
    print('Uploading series image...')
    set_file(api_url, api_key, series_image, response_ids, 'series_image')

    return 0

if __name__ == "__main__":
    main()
