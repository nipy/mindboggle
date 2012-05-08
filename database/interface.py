"""
Interface is used to establish a connection to MindBoggleDB and settings

"""
from bulbs.rest import Resource
from .base import Database


class Interface(object):
    """
    Establish the connection
    """
    def __init__(self,server):
        Database.resource = Resource(db_url=server)