from enum import IntEnum

class Axis(IntEnum):
    X = 0
    Y = 1

class Item:
    def __init__(self, id, dx, dy):
        self.Id = id
        self.Dx = dx
        self.Dy = dy
        self.Weight = dx * dy
        
    def __eq__(self, other):
        return self.Dx == other.Dx and self.Dy == other.Dy

    def __lt__(self, other):
        return ((self.Weight, self.Dx) < (other.Weight, other.Dx))

    def Dimension(self, axis):
        if axis == Axis.X:
            return self.Dx
        elif axis == Axis.Y:
            return self.Dy
        else:
            raise ValueError("Invalid axis.")

class Bin:
    def __init__(self, dx, dy):
        self.Dx = dx
        self.Dy = dy
        self.WeightLimit = dx * dy

    def Dimension(self, axis):
        if axis == Axis.X:
            return self.Dx
        elif axis == Axis.Y:
            return self.Dy
        else:
            raise ValueError("Invalid axis.")