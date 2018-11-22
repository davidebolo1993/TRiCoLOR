#!/usr/bin/python

#Plot precision and recall of a method with harmonic means (f1 scores)

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
		curve, = plt.plot(xs, ys, "--", color="gray", linewidth=0.5)
		plt.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="gray")



def plotPrecisionRecallDiagram(title=None, points=None, labels=None,colors=None):

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
		plt.legend(scps, labels, loc=0, scatterpoints=1, numpoints=1, fancybox=False) #loc=0 guess the best location for legend
	
  plt.axis([-0.02, 1.02, -0.02, 1.02])


#points must be a list of tuples [(pr1,re1),(pr2,re2)], labels a list-like ["A","B"], colors a list like ["#924242","#212977"]


plotPrecisionRecallDiagram("Title",points, labels, colors)
plt.show()


