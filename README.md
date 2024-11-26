# EMUCAT Components

The algorithms that a run in the [emucat workflows](https://github.com/ASKAP-EMUCat/emucat_workflows) pipeline. 

## Structure

emucat_components - science algorithms.
- component_1
    - component1_code_directory
    - Dockerfile
- component_2
    - component2_code_directory
    - Dockerfile
- component_n (not owned by project)
    - Dockerfile


Each component contains a `docker-compose.yml` produces a docker image which needs to pushed to the [AusSRC docker hub](https://hub.docker.com/orgs/aussrc/repositories) in order for the Nextflow pipelines to download them automatically.