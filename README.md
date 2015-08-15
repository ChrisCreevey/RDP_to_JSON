# RDP_to_JSON
tools for parsing RDP taxonomy output in JSON format for visualisation

This is a C programme to create a JSON file from a taxonomy assignment from the RDP classifier (https://rdp.cme.msu.edu/) combined with totals for abundances for the OTUs, in order to create a visualisations using the D3 javascript library (see https://github.com/mbostock/d3/wiki/Gallery for examples).

To compile the code simply run:

cc RDP_to_JSON.c -o RDP_to_JSON

(ignore the warnings)

Then run the tool using the provided examples with:

./RDP_to_JSON allrank_OTU.fas_classified.txt OTU_totals.txt 80

Where:

"allrank_OTU.fas_classified.txt" is the unchanged output from the allrank results from the RDP classifier

"OTU_totals.txt" is a file detailing the abundance of each OTU.

"80" represents a cutoff for the classification (representing 80%).

This produces a file called "taxonomy.json", which is used for visualisation.

The enclosed modified html files allow visualisation for the OTU data usign the provided html files.
