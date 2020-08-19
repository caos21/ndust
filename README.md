# ndust
Nanoparticle growth simulations in low-temperature plasmas.

- **crat**: precalculation for electrostatic and van der Waals nanoparticle interaction.
- **nevo**: evolution of nanoparticle growth. A straightforward python implementation is available [here](https://github.com/caos21/Grodi).
- **ndust**: Qt GUI.
- **press**: notebooks used for visualizations of the results.
- **test**: important unit tests.


## Requirements

### Ubuntu

- Install cmake, boost and sundials
 ```
 sudo apt install cmake g++ libboost-all-dev libsundials-serial-dev
 ```

## Compile
 
- Create build
```
cmake -Bbuild -H.
```
 
- Change to directory build
```
cd build
```
 
- Compile and link
```
make
