import sys

import networkx as nx
from networkx.algorithms.approximation import clique

from Model import Item, Bin
from PlacementPoints import PlacementPointStrategy

from SymmetryBreaking import SymmetryBreaking

class Preprocess:
    def __init__(self, items, bin, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        self.Items = items
        self.ProcessedItems = []
        self.Bin = bin

        self.IsCompleted = False

        self.UpperBoundsBin = sys.maxsize

        self.PlacementPointStrategy = placementPointStrategy
        #self.EnableNormalPatterns = True

        self.IncompatibleItems = set()
        self.RemovedItems = []
        self.FixItemToBin = []
        self.BinDomains = []

        self.ItemPlacementPatternsX = []
        self.ItemPlacementPatternsY = []
        self.GlobalPlacementPatternsX = []

    def DetermineConflicts(self, items, H, W):
        for i, itemI in enumerate(items):
            for j in range(i + 1, len(items)):
                itemJ = items[j]

                if itemI.Dx + itemJ.Dx > W and itemI.Dy + itemJ.Dy > H:
                    self.IncompatibleItems.add(frozenset((i, j)))
                    continue

    def RemoveLargeItems(self, items, H, W):
        filteredItemIndices = []
        for i, item in enumerate(items):
            dy = item.Dy
            dx = item.Dx

            if dy == H and dx == W:
                #print(f'Item {i} has the same dimensions as the bin and will be removed.')
                self.RemovedItems.append(Item(i, dx, dy))
                continue
            
            isFullyIncompatible = True
            for j, itemJ in enumerate(items):
                if i == j:
                    continue

                if item.Dx + itemJ.Dx > W and item.Dy + itemJ.Dy > H:
                    continue

                isFullyIncompatible = False
                break

            if isFullyIncompatible:
                #print(f'Item {i} is fully incompatible and will be removed.')
                self.RemovedItems.append(Item(i, dx, dy))
                continue

            filteredItemIndices.append(i)

        #self.BinPacking.Model.ObjCon = len(self.RemovedItems)

        newItems = []
        for index, i in enumerate(filteredItemIndices):
            newItems.append(Item(index, items[i].Dx, items[i].Dy))

        return newItems, len(filteredItemIndices)

    def Run(self):
        if self.IsCompleted:
            return 

        items = sorted(self.Items, reverse=True) # TODO: build conflict graph and compute maximal clique
        items, newNumberOfItems = self.RemoveLargeItems(items, self.Bin.Dy, self.Bin.Dx)

        self.DetermineConflicts(items, self.Bin.Dy, self.Bin.Dx)

        conflictGraph = nx.Graph()
        conflictGraph.add_nodes_from([item.Id for i, item in enumerate(items)])
        for i in range(len(items)):
            itemI = items[i]
            for j in range(len(items)):
                itemJ = items[j]
                if frozenset((itemI.Id, itemJ.Id)) in self.IncompatibleItems:
                    conflictGraph.add_edge(itemI.Id, itemJ.Id)

        maxClique = clique.max_clique(conflictGraph)
        sortedMaxClique = sorted(maxClique)

        newItems = []
        for i, oldIndex in enumerate(sortedMaxClique):
            newItems.append(Item(i, items[oldIndex].Dx, items[oldIndex].Dy))
        
        for i, item in enumerate(items):
            if item.Id in sortedMaxClique:
                continue
            newItems.append(Item(len(newItems), items[item.Id].Dx, items[item.Id].Dy))

        self.IncompatibleItems.clear()
        self.DetermineConflicts(newItems, self.Bin.Dy, self.Bin.Dx)

        newNumberOfItems = len(newItems)
        self.UpperBoundsBin = newNumberOfItems

        self.FixItemToBin = SymmetryBreaking.DetermineFixedItems(self.IncompatibleItems, newNumberOfItems)

        self.BinDomains = SymmetryBreaking.CreateReducedBinDomains(self.IncompatibleItems, newNumberOfItems, self.UpperBoundsBin, self.FixItemToBin)

        self.ItemPlacementPatternsX, self.ItemPlacementPatternsY, self.GlobalPlacementPatternsX = SymmetryBreaking.CreateBinDependentPlacementPatterns(self.IncompatibleItems, self.FixItemToBin, newItems, self.UpperBoundsBin, self.Bin, self.PlacementPointStrategy)

        self.ProcessedItems = newItems

        self.IsCompleted = True