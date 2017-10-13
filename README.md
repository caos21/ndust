# ndust
Nanoparticle growth simulations

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
