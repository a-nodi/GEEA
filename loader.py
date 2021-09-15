"""
A module that containing load functions

"""
import os
import json


def parse_json(json_dir: str) -> dict:
    """
    A function that reads json
    :param json_dir:
    :return:
    """
    with open(json_dir, 'r') as f:
        json_content = json.load(f)
        return json_content


def load_setting():
    """
    A function that reads setting.json and returns setting value dictionary
    :return dict_of_setting:
    """
    SETTING_JSON_DIR = os.path.join(os.getcwd(), 'resource\\json\\setting.json')
    dict_of_setting = parse_json(SETTING_JSON_DIR)
    return dict_of_setting


def load_resource_dir():
    """
    A function that reads resource_dir.json and returns resource_dir vaule dictionary
    :return dict_of_resource_dir:
    """
    RESOURCE_JSON_DIR = os.path.join(os.getcwd(), 'resource\\json\\resource_dir.json')
    dict_of_resource_dir = parse_json(RESOURCE_JSON_DIR)
    return dict_of_resource_dir
