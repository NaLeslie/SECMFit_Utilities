# SECMFit_Utilities
Utility methods for fitting arbitrary SECM features

## Notes about the structure:
- `run()` and `run2()` are dummy methods to stand-in for methods auto-generated by COMSOL Multiphysics. As is the `class Model`
- `comsolcompile` (the compiler that must be used for COMSOL api apps) does not like anything being outside of the main class that is not part of Java 1.7 jdk or their proprietary libraries. The `.java` file will compile, but the program always quietly crashes.
- This project is meant to just hold the code that is pasted into the auto-generated `.java` file from COMSOL. Some modifications need to be made within the `run()` and `run2()` methods. The contents of `/Example-impl` shows an example of this implementation.