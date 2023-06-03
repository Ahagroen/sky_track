# General
1. [ ] Refactor into modules for high level and low level logic
   1. [ ] Potentially, Track/Passes/Station as 3 crates all running sky_track for drivers
   2. [ ] Driver crate only for dependancies?
2. [ ] Ensure ease of integration into projects
   1. [ ] Build out API concept for each
   2. [ ] Multiple output forms, ensure serde support for the needed elements
3. [ ] Ensure that whole ecosystem interacts well together, including ground station drivers, websockets? etc.
4. [ ] Load/Allow for file based rulesets/preferences
   1. [ ] Station list
   2. [ ] Satellite lists
   3. [ ] Config

# Track
## Ground Track and Subpoints
1. [ ] Day/Night graphics
   1. [ ] MVP is satellite in/out of sunlight
   2. [ ] Support to compute in a general way for multiple rendering options
2. [ ] Support for JSON output
3. [ ] Refactor code into high level pub api and low level driver

# Passes
Segment handling satellite passes over a ground station, and associated ground segment propagation
## Pass List
1. [ ] Allow for multiple ground stations to be considered
   1. [ ] Loaded from YAML/.ini or other from disk - direct support
   2. [ ] Requested from whatever wrapper is calling it at runtime
2. [ ] Optimize pass finding algorithm
   1. [ ] Multithreading for multiple ground stations
3. [ ] Improve Output options
   1. [ ] JSON support

## Pass handling
1. [ ] Support for Generating pointing look angles 
   1. [ ] Real Time Update
   2. [ ] pass values cleanly for use by station
   3. [ ] Simple to use API wrapper
2. [ ] Consider handling station performance
   1. [ ] Data recieved
   2. [ ] Pointing error
   3. [ ] True pointing of the station
   4. [ ] Doppler
   5. [ ] Gain

