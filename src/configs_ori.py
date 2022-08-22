#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:47:55 2020

@author: tanvi
"""
import subprocess
import json
import os


class config:
    def __init__(self):
        self.runParams = {}
        self.conf = (
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            + "/configs/configs.json"
        )
        with open(self.conf) as file1:
            data = json.load(file1)
        for name in data["Machine"]:
            if name in subprocess.getoutput("hostname"):
                for key, val in data["Machine"][name].items():
                    self.runParams[key] = val
                self.runParams["host"] = "server"
            else:
                for key, val in data["local"].items():
                    self.runParams[key] = val
                self.runParams["host"] = "local"
