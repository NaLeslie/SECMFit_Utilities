# SECMFit_Utilities
Utility methods for fitting arbitrary SECM features

## Notes about the structure:
- `run()` and `runk()` are dummy methods to stand-in for methods auto-generated by COMSOL Multiphysics. As is the `Model` class.
- `comsolcompile` (the compiler that must be used for COMSOL api apps) does not like anything being outside of the main class that is not part of Java 1.7 jdk or their proprietary libraries. The `.java` file will compile, but the program always quietly crashes.
- `/src/` is meant to just hold the code that is copied into the auto-generated `.java` file from COMSOL. Some modifications need to be made within the `run()` and `runk()` methods. The contents of `/Example-impl` shows an example of this implementation.

## Operating systems:
- `/Example-impl/run.bat` will compile and run the simulation on Windows machines.
- On Windows systems, `C:\Program Files\COMSOL\COMSOL##\Multiphysics\bin\win64` must be added to the `PATH` environment variable.
- The comsol compiler has a different name in Unix-based systems (`comsole compile` instead of `comsolcompile`).
- COMSOL's security preferences need to be updated. In `File>Preferences>Security`:
  - Give COMSOL access to all files (so that you can write to files).
  - Give COMSOL access to system properties (so that you can find the current working directory).
