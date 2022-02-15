import os
import json

class BinSerializationNew:
    def __init__(self, dx, dy):
        self.Length = dx
        self.Height = dy

class ItemSerializationNew:
    def __init__(self, dx, dy, amount):
        self.Length = dx
        self.Height = dy
        self.Demand = amount

path = 'data/input/OPP/InterestingInstancesOld'
for root, dirs, files in os.walk(path):
    for fileName in files:
        if '.json' not in fileName:
            continue
        
        data = None
        filePath = os.path.join(path, fileName)
        with open(filePath, 'r') as file:
            data = file.read()
            data = json.loads(data)

        newJsonDict = {}
        newJsonDict["Name"] = data["Name"]

        binDx = data["Container"]["Length"]
        binDy = data["Container"]["Width"]
        bins = [BinSerializationNew(binDx, binDy)]

        items = []
        for item in data["ItemTypes"]:
            itemDx = item["Length"]
            itemDy = item["Width"]
            amount = item["Amount"]
            newItem = ItemSerializationNew(itemDx, itemDy, amount)

            items.append(newItem)

        newJsonDict["Objects"] = bins
        newJsonDict["Items"] = items

        with open(os.path.join("data/input/OPP/InterestingInstances", fileName), "w") as fp:
            print(f"Instance {fileName} saved")
            json.dump(newJsonDict, fp, indent=4, default=lambda x: x.__dict__)