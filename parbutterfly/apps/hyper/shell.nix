# nix expression to build and install Cilk+
let
  pkgs = import <nixpkgs> { };
  cilkgit = "https://gitlab.com/parallel-computing/cilkplus/cilkplus";
  cilkplus = pkgs.callPackage pkgs.example-cmake {
    name = "cilkplus";
    src = pkgs.fetchFromGitLab {
      url = cilkgit;
      rev = "2024-10-22";
      sha256 = "sha256-xxx"; # Replace with actual sha256 hash
    };
    buildInputs = [pkgs.cmake];
  };
in
  pkgs.stdenv.mkDerivation {
    name = "cilkplus-package";
    src = ./.;
    buildPhase = ''
      mkdir -p $out
      cp -r ${cilkplus} $out
    '';
    installPhase = ''
      install -DmR ${pkgs.cmake.dev} $out/
      install -DmR $out/ ${pkgs.cmake.dev}
      install -DmR ${pkgs.cmake.dev} $out/
      mkdir -p ${out}/include
      mkdir -p ${out}/lib
      cp -r $src/include ${out}/include
      cp -r $src/lib ${out}/lib
    '';
    outputs = "out";
  }