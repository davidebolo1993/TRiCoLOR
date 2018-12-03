#!/usr/bin/python

#Plot precision and recall with harmonic means (f1 scores) for couples of points

import scipy
import matplotlib.pyplot as plt

def fmeasure(p, r):

	"""Calculates the f1 score for precision p and recall r."""

	return float(2*p*r / (p+r))


def fmeasureCurve(f, p):

	"""For a given f1 value and precision get the recall value."""

	return float(f * p / (2 * p - f))


def plotFMeasures(fstepsize=.1, stepsize=0.001):

	"""Plots 10 f1 score curves."""

	p = scipy.arange(0., 1., stepsize)[1:]

	for f in scipy.arange(0., 1., fstepsize)[1:]:

		points = [(x, fmeasureCurve(f, x)) for x in p if 0 < fmeasureCurve(f, x) <= 1.5]
		xs, ys = zip(*points)
		curve, = plt.plot(xs, ys, "--", color="lightgray", linewidth=0.5)
		plt.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="lightgray")



def plotPrecisionRecallDiagram(title=None, points=None, labels=None,colors=None,loc=0):

	ax = plt.gca()   
	plt.title(title)
	plt.xlabel("Precision")
	plt.ylabel("Recall")
	plotFMeasures()
	scps = [] 
	
	for i, (x, y) in enumerate(points):
			
		label = labels[i]
		color=colors[i]
		scp = ax.scatter(x, y, label=label, c=color,s=50, linewidths=0.75, alpha=0.75)
		scps.append(scp)
		plt.legend(scps, labels, loc=loc, scatterpoints=1, numpoints=1, fancybox=False) #loc=0 guess the best location for legend
	
	plt.axis([-0.02, 1.02, -0.02, 1.02])


#points must be a list of tuples, such as: [(pr1,re1),(pr2,re2),...];
#labels are a list-like object, such as: ["A","B",...];
#colors are a list-likeobject, such as: ["#924242","#212977",...]
#loc must be an integer or a tuple of coordinates, as specified here: https://matplotlib.org/api/legend_api.html


plotPrecisionRecallDiagram("Title",points, labels, colors)
plt.show()



# Plot precision and recall with harmonic means (f1 scores) for groups of cuples of points: points of each group are connected by lines


import scipy
import matplotlib.pyplot as plt

def fmeasure(p, r):

	"""Calculates the f1 score for precision p and recall r."""

	return float(2*p*r / (p+r))


def fmeasureCurve(f, p):

	"""For a given f1 value and precision get the recall value."""

	return float(f * p / (2 * p - f))


def plotFMeasures(fstepsize=.1, stepsize=0.001):

	"""Plots 10 f1 score curves."""

	p = scipy.arange(0., 1., stepsize)[1:]

	for f in scipy.arange(0., 1., fstepsize)[1:]:

		points = [(x, fmeasureCurve(f, x)) for x in p if 0 < fmeasureCurve(f, x) <= 1.5]
		xs, ys = zip(*points)
		curve, = plt.plot(xs, ys, "--", color="lightgray", linewidth=0.5)
		plt.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="lightgray")



def plotPrecisionRecallDiagram(title=None, points=None, labels=None,colors=None,linecolors=None,loc=0):

	ax = plt.gca()   
	plt.title(title)
	plt.xlabel("Precision")
	plt.ylabel("Recall")
	plotFMeasures()
	scps = [] 

	for h in range(len(points)):

		for i, (x,y) in enumerate(points[h]):
			
			label = labels[h][i]
			color=colors[h][i]
			scp = ax.scatter(x, y, label=label, c=color,s=50, linewidths=0.75, alpha=0.75)
			scps.append(scp)

		if linecolors is not None:

			lines = ax.plot([el[0] for el in points[h]],[el[1] for el in points[h]],":", color=linecolors[h], linewidth=0.5)

		else:

			lines= ax.plot([el[0] for el in points[h]],[el[1] for el in points[h]],":", color="black", linewidth=0.5)
		
		plt.legend(loc=loc, scatterpoints=1, numpoints=1, fancybox=False) #loc=0 guess the best location for legend
	
	plt.axis([-0.02, 1.02, -0.02, 1.02])


#points must be a list of lists of tuples, such as: [[(pr1.0,re1.0),(pr2.0,re2.0),...],[(pr1.1,re1.1),(pr2.1,re2.1),...],...]
#labels are a list of list-like object, such as: [["A","B",...],["C","D",...],...]blu
#colors are a list of list-like object, such as: [["#924242","#212977",...],["#23D23F", "#EEAA12",...],...]
#example of blue and red color palettes with 5 points that can be used:
#blue color palette:["#ece7f2", "#a6bddb","#2b8cbe", "#4522D8", "#1C1960"]
#red color palette:["#fee8c8", "#fdbb84", "#e34a33", "#E21919", "#B01212"]
#loc must be an integer or a tuple of coordinates, as specified here: https://matplotlib.org/api/legend_api.html


plotPrecisionRecallDiagram("Precision and Recall for exact number of repetitions",result, labels, colors, linecolors=["blue", "red"])
plt.show()
