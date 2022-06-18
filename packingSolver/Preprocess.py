from re import A
import sys
from typing import ItemsView
import math

import networkx as nx
from networkx.algorithms.approximation import clique

from Model import Item, Bin

from OrthogonalPacking import OrthogonalPackingSolver
from SymmetryBreaking import SymmetryBreaking
#from PlacementPoints import PlacementPointStrategy
import PlacementPoints

""" 
Implements the reduction procedures from section 3 in 
Clautiaux, F., Carlier, J., & Moukrim, A. (2007). A new exact method for the two-dimensional orthogonal packing problem. 
European Journal of Operational Research, 183(3), 1196-1211. 
"""
class PreprocessOrthogonalPacking:
    def __init__(self, items, bin):
        self.Items = items
        self.Bin = bin
        self.PreprocessedItems = items
        self.PreprocessBin = bin

    def Run(self):
        newItems = list(self.Items)
        newBin = Bin(self.Bin.Dx, self.Bin.Dy)

        self.RemoveLargeItems(newItems, newBin)
        self.FilterFrameConfiguration(newItems, newBin)
        self.EnlargeItemsByOrthogonalPacking(newItems, newBin)

        if len(self.Items) > len(newItems):
            print(f"Preprocess removed {len(self.Items) - len(newItems)} items and reduced container area by {self.Bin.Dx * self.Bin.Dy - newBin.Dx * newBin.Dy}.")

        self.PreprocessedItems = newItems
        self.PreprocessBin = newBin

    def RemoveLargeItems(self, items, bin):
        itemRemoved = True
        while itemRemoved:
            itemRemoved = False
            itemsToRemove = []
            for i, item in enumerate(items):
                if item.Dx == bin.Dx:
                    bin.Dy -= item.Dy

                    itemsToRemove.append(item)
                    itemRemoved = True

                    break

                if item.Dy == bin.Dy:
                    bin.Dx -= item.Dx

                    itemsToRemove.append(item)
                    itemRemoved = True

                    break

            for item in itemsToRemove:
                items.remove(item)

    def FilterFrameConfiguration(self, items, bin):
        for i, itemI in enumerate(items):
            for j, itemJ in enumerate(items):
                if j == i:
                    continue

                for k, itemK in enumerate(items):
                    if k == j or k == i:
                        continue

                    for l, itemL in enumerate(items):
                        if l == k or l == j or l == i:
                            continue

                        if (itemI.Dx + itemL.Dx == itemJ.Dx + itemK.Dx and
                            itemI.Dx + itemL.Dx == bin.Dx and
                            itemI.Dy + itemK.Dy == itemJ.Dy + itemL.Dy and
                            itemI.Dy + itemK.Dy == bin.Dy):
                            bin.Dx = bin.Dx - itemI.Dx - itemJ.Dx
                            bin.Dy = bin.Dy - itemK.Dy - itemL.Dy

                            itemsToRemove = set(itemI, itemJ, itemK, itemL)

                            for item in itemsToRemove:
                                items.remove(item)

                            return

    def EnlargeItemsByOrthogonalPacking(self, items, bin):
        feasibleThresholdToItemCount = {}
        largestItemThreshold = -1
        largestItemCount = 0

        for p in range(1, math.floor(bin.Dy / 2.0) + 1):
            
            largeItemSubset = []
            smallItemSubset = []
            itemSubset = []

            largeItemDx = 0
            largestSmallItemDx = 0
            
            for i, item in enumerate(items):
                if item.Dy >= bin.Dy - p:
                    largeItemDx += item.Dx

                    itemSubset.append(item)
                    largeItemSubset.append(item)
                    continue

                if item.Dy <= p:
                    largestSmallItemDx = max(largestSmallItemDx, item.Dx)

                    itemSubset.append(item)
                    smallItemSubset.append(item)
                    continue
            
            if len(largeItemSubset) == 0 or len(smallItemSubset) == 0 or len(itemSubset) == 0 or largestSmallItemDx > largeItemDx:
                continue

            reducedBin = Bin(largeItemDx, bin.Dy)
            
            print(f"Preprocess OPP call with {len(itemSubset)} / {len(items)}")

            orthogonalPackingSolver = OrthogonalPackingSolver(itemSubset, reducedBin)
            isFeasible = orthogonalPackingSolver.Solve(False)
            
            if isFeasible:
                if len(itemSubset) > largestItemCount:
                    largestItemCount = len(itemSubset)
                    largestItemThreshold = p

                    feasibleThresholdToItemCount[p] = itemSubset

                    if largestItemCount == len(items):
                        # All items could be removed.
                        break

        if largestItemThreshold == -1:
            return

        maximalItemSet = feasibleThresholdToItemCount[largestItemThreshold]
        itemsToRemove = []

        reducedBinDx = 0
        for item in maximalItemSet:
            itemsToRemove.append(item)
            reducedBinDx += item.Dx

        bin.Dx -= reducedBinDx

        for item in itemsToRemove:
            items.remove(item)

class PreprocessBinPacking:
    def __init__(self, items, bin, placementPointStrategy = PlacementPoints.PlacementPointStrategy.UnitDiscretization):
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