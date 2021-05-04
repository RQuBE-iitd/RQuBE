# RQuBE

## Introduction
Regular path query (RPQ) on labelled graphs is one of the core operators in graph analytics. In the decision version of an RPQ,the input is a pair of source-destination nodes, a labelled graph, and a regular expression. The goal is to determine if there exists a simple path from the source to destination, such that the label sequence in the path satisfies the input regular expression. However, identifying this regular expression is a tedious task and the goal of our work isto remove the barrier imposed by regular expressions without compromising on their expressive power. 
Thus We use query by example paradigm where the user provides an exemplar source-destination pair. The exemplar pair acts as a proxy for the regex constraint and communicates to the query evaluation engine the constraints in a more user-friendly manner. However,to make sense of this query, the query evaluation engine needs to execute the following tasks:
(1) Infer the regular expression that characterizes the paths between the source and the destination, and
(2) Identify nodes in the graph that are connected in a similar manner to the source as characterized by the regex inferred from the exemplar pair.

## Datasets
Datasets can be downloaded from here: https://drive.google.com/drive/folders/10odFg2fhUzrXc6symOcWQf9-abcW1n3k?usp=sharing

## Running experiments
To run a dataset with edge labels use x_edge file otherwise x_node file.

### BBFS
