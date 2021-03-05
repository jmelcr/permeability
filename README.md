# Permeability of small polar solutes through lipid membranes at various phases

This is a project that started independently 
as my work with permeability of ethanol trhough lipid membranes and 
the work in Poolman's lab done by Jacopo Frallicciardi. 

They have developed an accurate technique to measure permeability of small polar compounds 
and now they are interested in understanding how such compounds, especially weak acids and water,
permeate through lipid membranes in the context of yeast lipids and yeast plasma membrane. 

My simulations are instrumental in understanding what are the driving forces in the process and
what happens on the molecular level. 

Long ago since last time!

# Current status

We have advanced with the msp with Jacopo. 
The (dis)agreement between sim/exp lead us to adopt the "Small" size 
as the best candidate for the permeating particle. 

Now, all results are generated using this particle size. 

The results include series of:
 - varying tail chain lengths
 - varying sterol concentrations 
 - various saturation levels (PO, DP, + intermed.)
 - newly, we have also decied to include permants of varying hydrophobicity ...

## particles with various hydrophilicity/phobicity

After the first read by Bert, 
he chose to include back a figure comparing the hydrophilicity of compounds (LogPOW) with their permeability (logPerm).
This is called the Meyer-Overton rule and hold well for lots of rather hydrophobic compounds. 

In our work, we find that the hydrophilic compounds (weak acids that partition better to water that to oils and e.g. glycerol)
have a different, steeper, slope of the dependence of permeability on their particion coefficient. 

To probe this effect, I have used a series of small particles from SP3 to SC5.
The different slope is recovered from simulations 
supporting the idea that solubility is the key feature that determines the permeability. 
The difference in the slope seems to arise mainly from the capacity of the hydrophilic compounds
to attract water molecules deeper inside the membrane core,
which costs extra energy. 
The more hydrophilic, the more water dragged inside and the more the cost for permeation. 

This part of research is contained in the latest commits above or below this message. 
Its placement in in the folder structure is not the most optimal, but let's work with it as an exception for the time being. 
This is/will be noted in the corresponding Readme file. 
