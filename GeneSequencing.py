#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4':
# 	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH_COST = -3
INDEL_COST = 5
SUB_COST = 1
MATCH_OPERATION = "Match"
INDEL_OPERATION = "Insertion/Deletion"
SUB_OPERATION = "Substitution"
D = 3
TOP = "Top"
LEFT = "Left"
DIAGONAL = "Diagonal"
INFINITY = float('inf')

class GeneSequencing:

	def __init__( self ):
		pass

	def buildUnrestrictedTable(self, numRows, numCols):
		thing = ((None, None), (None, None))  # ((cost, operation),(parentRow, parentCol))

		if numRows > self.MaxCharactersToAlign:
			numRows = self.MaxCharactersToAlign
		if numCols > self.MaxCharactersToAlign:
			numCols = self.MaxCharactersToAlign

		self.table = [[thing] * (numCols + 1) for g in range(numRows + 1)]

	def buildBandedTable(self, numRows, numCols):

		# Build a k by n array of arrays
		thing = ((None, None), (None, None))  # ((cost, operation),(parentRow, parentCol))
		if numRows > self.MaxCharactersToAlign:
			numRows = self.MaxCharactersToAlign
		if numCols > self.K:
			numColsforInit = self.K
		table = [[thing] * numColsforInit for g in range(numRows + 1)]

		# Delete the elements that we don't need from the front rows in the table
		itemsToDelete = D
		lastTableRowIndex = len(table) - 1
		for i in range(D + 1):
			for j in range(itemsToDelete):
				del table[i][0]
				del table[lastTableRowIndex - i][0]
			itemsToDelete -= 1

		self.table = table

		# breakingPointColIndex= numCols - D
		# breakingPointRowIndex = numRows - D
		# table = [[None]]
		# diagonalAnchor = 0
		# done = False
		#
		# d = D
		#
		# while not done:
		# 	# Add d values to your right
		# 	for i in range(d):
		# 		# Account for base cases while starting the table
		# 		if (len(table) == 1):
		# 			table[diagonalAnchor].append(i * INDEL_COST)
		# 		else:
		# 			table[diagonalAnchor].append((None))
		# 	# Add D values below you
		# 	if diagonalAnchor == 0:
		# 		# We need to create the next D rows so we can add to them
		# 		for i in range(d):
		# 			table.append(i * INDEL_COST) # account for the base case
		# 	elif diagonalAnchor > 0 and diagonalAnchor < (breakingPointRowIndex + 1):
		# 		# The next D - 1 rows have already been added
		# 		# Add just one row
		# 		table.append([None])
		# 	# else:
		# 	# Your rows have already been added. So do nothing.
		#
		# 	indexOfLastRowInTable = len(table) - 1
		# 	indexOfLastColInTable = len(table[indexOfLastRowInTable]) - 1
		#
		# 	if indexOfLastRowInTable >= breakingPointRowIndex and indexOfLastColInTable >= breakingPointColIndex:
		# 		d -= 1
		#
		# self.table = table

	def setBestCost(self, row, col):
		if self.banded:
			topCost, topPos = self.getTopValAndPos(self.table, row, col)
			leftCost, leftPos = self.getLeftAndPos(self.table, row, col)
			diagonalCost, diagPos = self.getDiagonalValAndPos(self.table, row, col)
		else:
			# Add the base case for the first row
			if col == 0 and row == 0:
				self.table[row][col] = ((0, None), (None, None))
				return
			elif col == 0:
				self.table[row][col] = ((INDEL_COST * row, INDEL_OPERATION), (None, None))
				return
			elif row == 0:
				self.table[row][col] = ((INDEL_COST * col, INDEL_OPERATION), (None, None))
				return

			leftCost = self.table[row][col - 1][0][0]
			leftPos = (row, col - 1)
			topCost = self.table[row - 1][col][0][0]
			topPos = (row - 1, col)
			diagonalCost = self.table[row - 1][col - 1][0][0]
			diagPos = (row - 1, col - 1)

		# If the two chars match in the gene sequence...
		if self.areCharsEqual(row, col):
			cost = diagonalCost + MATCH_COST
			operation = MATCH_OPERATION
			if self.banded:
				junk, parentLocation = self.getDiagonalValAndPos(self.table, row, col)
			else:
				parentLocation = (row - 1, col - 1)
		else:
			if leftCost <= topCost and leftCost <= diagonalCost: # If leftCost is the smallest
				cost = leftCost + INDEL_COST
				operation = INDEL_OPERATION
				parentLocation = leftPos
			elif topCost < leftCost and topCost <= diagonalCost: # If topCost is the smallest
				cost = topCost + INDEL_COST
				operation = INDEL_OPERATION
				parentLocation = topPos
			elif diagonalCost < leftCost and diagonalCost < topCost: # If diagonalCost is the smallest
				cost = diagonalCost + SUB_COST
				operation = SUB_OPERATION
				parentLocation = diagPos
			else:
				print("ERROR. You should never get here.")
		self.table[row][col] = ((cost, operation), parentLocation)
		# return table

	def areCharsEqual(self, row, col):
		posInSeq1 = row - 2 # Nothing is tricky with the row index
		isBanded = self.banded
		if not isBanded or row <= D:
			posInSeq2 = col - 2
		else:
			actualPos = row - D - 2
			posInSeq2 = actualPos

		# Compare string characters
		if self.seq1[posInSeq1] == self.seq2[posInSeq2]:
			return True
		return False

	def getTopValAndPos(self, row, col):
		isLastColumn = False
		if (len(self.table[row]) - 1) == col:
			isLastColumn = True
		if len(row) <= self.K and row <= D and not isLastColumn:
			top = self.table[row - 1][col][0][0]
			pos = (row - 1, col)
		elif len(row) == self.K and row > D and row <= (len(self.table) - (D + 1)) and not isLastColumn:
			top = self.table[row - 1][col + 1][0][0]
			pos = (row - 1, col + 1)
		elif row > len(self.table) - (D + 1):
			top = self.table[row - 1][col + 1][0][0]
			pos = (row - 1, col + 1)
		else:
			top = INFINITY
			pos = (None, None)
		return top, pos

	def getLeftValAndPos(self, row, col):
		if col == 0: # The item has no left neighbor
			return INFINITY, (None, None)
		else:
			left = self.table[row][col - 1][0][0]
			pos = (row, col - 1)
			return left, pos

	def getDiagonalValAndPos(self, row, col):
		if row <= D:
			diagonal = self.table[row - 1][col - 1][0][0]
			pos = (row - 1, col - 1)
		else:
			diagonal = self.table[row - 1][col][0][0]
			pos = (row - 1, col)
		return diagonal, pos

	def backTraverse(self, row, col, returnString):
		if row == 0 and col == 0:
			return -1, -1, returnString
		operation = self.table[row][col][0][1]
		if operation == INDEL_OPERATION:
			returnString = '-' + returnString
			leftVal, leftPos = self.getLeftValAndPos(row, col)
			if self.banded:
				topVal, topPos = self.getTopValAndPos(row, col)
			else:
				topVal = self.table[row - 1][col][0][0]
				topPos = (row - 1, col)
			if topVal < leftVal:
				row = topPos[0] # Accessing a tuple
				col = topPos[1]
			else:
				row = leftPos[0] # Accessing a tuple
				col = leftPos[1]
		elif operation == MATCH_OPERATION or operation == SUB_OPERATION:
			if self.banded:
				diagVal, diagPos = self.getDiagonalValAndPos(row, col)
				row = diagPos[0]
				col = diagPos[1]
			else:
				row -= 1
				col -= 1
			lastPosInSeq1 = len(self.seq1) - 1
			character = self.seq1[lastPosInSeq1]
			self.seq1 = self.seq1[:lastPosInSeq1] # delete last char in seq1
			returnString = character + returnString
		return row, col, returnString

	def getNextColumnStartingPosition(self, row):
		if row >= D:
			return 0
		return 1

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length): # align_length. Number of characters they want us return in the align gene sequence.
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.seq1 = seq1
		self.seq2 = seq2
		self.K = (2 * D) + 1
		numRows = len(self.seq1)
		numCols = len(self.seq2)

		# Build the table before traversing it
		if(banded):
			self.buildBandedTable(numRows, numCols)
			for rowIndex in range(len(self.table)):
				for colIndex in range()
		else:
			self.buildUnrestrictedTable(numRows, numCols)
			for rowIndex in range(len(self.table)):
				for colIndex in range(len(self.table[rowIndex])):
					self.setBestCost(rowIndex, colIndex)

		# Table has been populated. Time for back traversal
		row = len(self.table) - 1
		col = len(self.table[row]) - 1
		returnString = ""

		while(row > -1 and col > -1):
			row, col, returnString = self.backTraverse(row, col, returnString)

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		lastRowIndex = len(self.table) - 1
		lastColIndex = len(self.table[lastRowIndex]) - 1
		score = self.table[lastRowIndex][lastColIndex][0][0]
		alignment1 = returnString[:align_length] + '  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = returnString[:align_length] + '  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


