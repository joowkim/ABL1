#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/1/22

@author: jaykim
"""

import json
import os
import subprocess
from typing import Dict


class Config:
    def __init__(self):
        self.run_params: Dict = dict()
        self.conf: str = (
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            + "/configs/configs.json"
        )
        self.hostname: str = subprocess.getoutput("hostname").strip()

        assert os.path.isfile(self.conf), f"{self.conf} is not found!"

        with open(self.conf) as fin:
            data = json.load(fin)

        for name in data["Machine"]:
            # this is for getting the server info
            if name in self.hostname:
                for key, val in data["Machine"][name].items():
                    self.run_params[key] = val
                self.run_params["host"] = "server"

            # this is for getting the local pc info
            else:
                for key, val in data["local"].items():
                    self.run_params[key] = val
                self.run_params["host"] = "local"
