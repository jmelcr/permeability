created from the simulations with tiny particle by copying only topologies 
and representative configurations and awh-profiles for gel/fluid phases (manually)

Only changed the particle in toplogy by simply

    find -name *.top | xargs sed 's/EOLT/EOLR/g' -i 
