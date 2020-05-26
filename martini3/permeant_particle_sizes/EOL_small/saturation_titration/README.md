created from the simulations with regular (R) particle by simply copying what was done there before running any simulations
including representative configurations and awh-profiles for gel/fluid phases (manually)

Only changed the particle in toplogy by simply

    find -name *.top | xargs sed 's/EOLR/EOLS/g' -i 
