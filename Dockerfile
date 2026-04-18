# syntax=docker/dockerfile:1.7
#
# ferx: R package + Rust (ferx-nlme) engine + RStudio Server
#
# Build (from this directory):
#   docker build -t ferx:latest .
#
# Run RStudio Server:
#   docker run --rm -p 8787:8787 -e PASSWORD=ferx ferx:latest
#   -> http://localhost:8787   user: rstudio   password: ferx
#
# Run the ferx CLI directly:
#   docker run --rm -v "$PWD:/work" -w /work ferx:latest ferx model.ferx --data data.csv

FROM rocker/tidyverse:latest

# ---------------------------------------------------------------------------
# 1. System build + runtime deps. Kept slim; LLVM/clang are installed later
#    only for the Enzyme plugin build and purged afterwards.
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake ninja-build lld libssl-dev pkg-config \
        python3 build-essential curl git ca-certificates \
        libzstd-dev wget gnupg lsb-release \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. rustup + upstream nightly. CARGO_HOME/RUSTUP_HOME under /opt so the
#    `rstudio` user (default RStudio Server user) can use the toolchain too.
# ---------------------------------------------------------------------------
ENV CARGO_HOME=/opt/cargo \
    RUSTUP_HOME=/opt/rustup \
    PATH=/opt/cargo/bin:/usr/local/bin:/usr/bin:/bin

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain nightly --profile minimal \
    && rustup component add rust-src --toolchain nightly \
    && chmod -R a+rX /opt/cargo /opt/rustup

# ---------------------------------------------------------------------------
# 3. Build the Enzyme LLVM plugin against nightly's matching LLVM, drop it
#    into the rustc sysroot, and register as the `enzyme` toolchain.
#    Everything in one RUN layer so the ~5 GB build tree never persists.
#    Follows the 30-min shortcut from INSTALL-SYSADMIN.md.
# ---------------------------------------------------------------------------
RUN set -eux; \
    LLVM_MAJOR=$(rustc +nightly --version --verbose \
        | sed -n 's/^LLVM version: \([0-9]*\).*/\1/p'); \
    wget -qO /tmp/llvm.sh https://apt.llvm.org/llvm.sh; \
    chmod +x /tmp/llvm.sh; \
    /tmp/llvm.sh "${LLVM_MAJOR}"; \
    chmod 644 /etc/apt/trusted.gpg.d/apt.llvm.org.asc || true; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        "llvm-${LLVM_MAJOR}-dev" "clang-${LLVM_MAJOR}"; \
    \
    git clone --depth 1 https://github.com/EnzymeAD/Enzyme /tmp/enzyme-src; \
    cmake -G Ninja -S /tmp/enzyme-src/enzyme -B /tmp/enzyme-build \
        -DLLVM_DIR="/usr/lib/llvm-${LLVM_MAJOR}/lib/cmake/llvm" \
        -DENZYME_CLANG=OFF -DENZYME_FLANG=OFF; \
    cmake --build /tmp/enzyme-build; \
    \
    TARGET=$(rustc +nightly -vV | sed -n 's/^host: //p'); \
    SYSROOT=$(rustc +nightly --print sysroot); \
    install -m 0644 \
        "/tmp/enzyme-build/Enzyme/LLVMEnzyme-${LLVM_MAJOR}.so" \
        "${SYSROOT}/lib/rustlib/${TARGET}/lib/libEnzyme-${LLVM_MAJOR}.so"; \
    \
    rustup toolchain link enzyme "${SYSROOT}"; \
    rustup default enzyme; \
    \
    rm -rf /tmp/enzyme-src /tmp/enzyme-build /tmp/llvm.sh; \
    apt-get purge -y --auto-remove "llvm-${LLVM_MAJOR}-dev" "clang-${LLVM_MAJOR}"; \
    rm -rf /var/lib/apt/lists/*

ENV RUSTUP_TOOLCHAIN=enzyme \
    RUSTFLAGS="-Z autodiff=Enable"

# ---------------------------------------------------------------------------
# 4. Clone ferx-nlme from GitHub, build the `ferx` CLI binary, keep the
#    source tree at /opt/ferx-nlme (the R package's Cargo.toml has a
#    relative path dep: `path = "../../../ferx-nlme"`). Drop build artifacts
#    and cargo caches so the layer stays small.
# ---------------------------------------------------------------------------
RUN set -eux; \
    git clone --depth 1 https://github.com/InsightRX/ferx-nlme /opt/ferx-nlme; \
    cd /opt/ferx-nlme; \
    cargo build --release; \
    install -m 0755 target/release/ferx /usr/local/bin/ferx; \
    rm -rf target .git; \
    rm -rf /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git

# ---------------------------------------------------------------------------
# 5. Copy the R package source, install it (which rebuilds the Rust staticlib
#    against the retained /opt/ferx-nlme via the relative path dep in
#    src/rust/Cargo.toml), then clean all build/caching state.
#    `remotes` is already in the rocker/tidyverse base image.
# ---------------------------------------------------------------------------
COPY . /opt/ferx
RUN set -eux; \
    R -e "remotes::install_local('/opt/ferx', dependencies=TRUE, upgrade='never')"; \
    rm -rf /opt/ferx/src/rust/target \
           /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git \
           /tmp/Rtmp* /tmp/downloaded_packages

# ---------------------------------------------------------------------------
# 6. Make the toolchain visible to interactive R sessions inside RStudio,
#    so `remotes::install_github("InsightRX/ferx")` etc. work from the UI.
# ---------------------------------------------------------------------------
RUN printf 'PATH=/opt/cargo/bin:${PATH}\nRUSTUP_TOOLCHAIN=enzyme\n' \
        > /home/rstudio/.Renviron \
    && chown rstudio:rstudio /home/rstudio/.Renviron

# RStudio Server entrypoint, port 8787, default user `rstudio` are all
# inherited from rocker/tidyverse — nothing more to do here.
