# Web API Interface

Base url : /API

## Input part related services

*/uploadGeolInput/*
POST .input file
Request format :
    Content-Type : multipart/form-data,
    Name of the part in the body : file
Let the user upload a geologic.input file, needed for the simulation.

*/uploadHydroInput/*
POST .input file
Request format :
    Content-Type : multipart/form-data,
    Name of the part in the body : file
Let the user upload a hydrologic.input file, needed for the simulation.

*/uploadMorphoInput/*
POST .input file
Request format :
    Content-Type : multipart/form-data,
    Name of the part in the body : file
Let the user upload a morphologic.input file, needed for the simulation.

*/uploadParamInput/*
POST .param file
Request format :
    Content-Type : multipart/form-data,
    Name of the part in the body : file
Let the user upload the param file, needed for the simulation.

## Simulation part related services
*/simulate/*
GET with no param, return nothing. Uses the HTTP response code to indicate
whether the simulation run well or not.
Launch the simulation

*/missingElementToSimulate/*
GET, return a JSON array
Give a list of the missing elements needed in order to run the simulation.

## Simulation results part related services
*/results/*
GET, return a JSON array
Give the list of results from the simulation.

*/resultsFormat/*
GET, return a JSON object
The list of the variables that will be produced as the result of the simulation.
