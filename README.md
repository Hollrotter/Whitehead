# Whitehead

The library Whitehead¹ is developted for the efficient aeroelastic analysis of Membrane Structures.

## Features

### Structural Analysis of Membranes

Implemented and validated. The theory with the verification and validation is published in a [Journal Paper](https://link.springer.com/article/10.1007/s00466-026-02810-w).

### Aerodynamic Analysis of thin surfaces based on Potential Theory

Under development

### Coupling of Structure and Aerodynamics

Not yet implemented

### Beams in contact with the Membrane Structure

Not yet implemented

## Building

### 1. Clone the repository

```
$ git clone --recursive https://github.com/Hollrotter/Whitehead
$ cd Whitehead
```

### 2. Compile the library

```
$ make
```

### 3. Run one of the test files in VSCode (tested with Linux)

```
$ code .
```

The following are example files for running a test file in VSCode:

c_cpp_properties.json
```
{
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
                "${workspaceFolder}/**",
                "${workspaceFolder}/../include"
            ],
            "defines": [],
            "compilerPath": "/usr/bin/gcc",
            "intelliSenseMode": "linux-gcc-x64"
        }
    ],
    "version": 4
}
```

tasks.json:
```
{
	"version": "2.0.0",
	"tasks": [
		{
			"type": "cppbuild",
			"label": "C/C++: g++ build active file",
			"command": "/usr/bin/g++",
			"args": [
				"-fdiagnostics-color=always",
				"-g",
				"-Ofast",
				"-Wall",
				"${file}",
				"-I${fileDirname}/../include",
				"-L${fileDirname}/../lib",
				"-l:libWhitehead.a",
				"-larmadillo",
				"-fopenmp",
				"-std=c++23",
				"-o",
				"${fileDirname}/test"
			],
			"options": {
				"cwd": "${fileDirname}"
			},
			"problemMatcher": [
				"$gcc"
			],
			"group": {
				"kind": "build",
				"isDefault": true
			},
			"detail": "compiler: /usr/bin/g++"
		}
	]
}
```

launch.json
```
{
    "configurations": [
        {
            "name": "C/C++: g++ build and debug active file",
            "type": "cppdbg",
            "request": "launch",
            "program": "${fileDirname}/test",
            "args": [],
            "stopAtEntry": false,
            "cwd": "${fileDirname}",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                },
                {
                    "description": "Set Disassembly Flavor to Intel",
                    "text": "-gdb-set disassembly-flavor intel",
                    "ignoreFailures": true
                }
            ],
            "preLaunchTask": "C/C++: g++ build active file",
            "miDebuggerPath": "/usr/bin/gdb"
        }
    ],
    "version": "2.0.0"
}
```

## Build Dependencies

* C++ compiler (g++, C++23)
* openMP
* Armadillo

[1]: The library is named after Gustave Whitehead (1874 - 1927) who was a pioneer of early aviation. Whitehead's Nr. 21 will be a key validation example.