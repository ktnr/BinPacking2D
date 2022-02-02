import json
import math

import os, sys
import re

import Model

import numpy as np

class BinSerialization:
    def __init__(self, dx, dy, dz, capacity = -1):
        self.Capacity = capacity
        self.Length = dx
        self.Width = dy
        self.Height = dz

class ItemSerialization:
    def __init__(self, id, dx, dy, dz, quantity):
        self.Id = id
        self.Volume = dx * dy * dz
        self.Length = dx
        self.Width = dy
        self.Height = dz
        self.Amount = quantity

class Converter:
    @staticmethod
    def ConvertFromSubproblem(items, container, timeLimit, instanceId):
        newContainer = BinSerialization(container.Dx, container.Dy, 1)
        newItems = []

        areaSum = 0.0
        for i, item in enumerate(items):
            newItem = ItemSerialization(i, item.Dx, item.Dy, 1, 1)

            areaSum += item.Dx * item.Dy

            newItems.append(newItem)

        numberOfItemTypes = len(newItems)

        epsilon = 100 * (1.0 - float(areaSum) / float(container.Dx * container.Dy))

        itemWidths = np.array([item.Dx for item in items])
        itemHeights = np.array([item.Dy for item in items])

        percentileDx25 = int(np.percentile(itemWidths, 25))
        percentileDx50 = int(np.percentile(itemWidths, 50))
        maxDx = max(itemWidths)

        percentileDy25 = int(np.percentile(itemHeights, 25))
        percentileDy50 = int(np.percentile(itemHeights, 50))
        maxDy = max(itemHeights)

        newJsonDict = {}

        name = f"{instanceId}e{int(epsilon)}w{container.Dx}h{container.Dy}n{len(items)}dx-25-50-100_{percentileDx25}_{percentileDx50}_{maxDx}dy-25-50-100_{percentileDy25}_{percentileDy50}_{maxDy}t{int(timeLimit)}"
        newJsonDict["Name"] = name
        newJsonDict["InstanceType"] = "2D-OPP"
        newJsonDict["NumberItemTypes"] = numberOfItemTypes

        newJsonDict["Container"] = newContainer
        newJsonDict["ItemTypes"] = newItems

        with open(os.path.join("data", "2D-OPP-subproblems", name + ".json"), "w") as fp:
            print(f"Instance {name} saved")
            json.dump(newJsonDict, fp, indent=4, default=lambda x: x.__dict__)