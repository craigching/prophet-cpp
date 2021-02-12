## Pre-requisites

```
$ brew install catch2
```

## Building

Make sure and clone with submodules

```
git clone --recursive https://github.com/craigching/prophet-cpp.git
```

Build `cmdstan`

```
$ cd cmdstan
$ make stan-update
$ make build
$ cd ..
```

Build prophet

```
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
$ cd ..
```
