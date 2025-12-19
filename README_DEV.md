# Developer Guide: Modern Python Project Template

[阅读简体中文版](i18n/README_zh.md)

This document is for developers using this project template. It explains the philosophy, tools, and workflows for building a high-quality Python application.

## Core Tech Stack

This template is built upon a foundation of modern, high-performance tools:

- **`uv`**: An extremely fast Python package installer and resolver, used for all dependency management, environment creation, and task execution.
- **`Ruff`**: An extremely fast Python linter and code formatter.
- **`Mypy`**: The standard for static type checking.
- **`Pytest`**: The framework for robust and scalable testing.
- **`Hatch`**: The build backend, with versioning managed by `hatch-vcs`.
- **GitHub Actions**: For automated CI/CD, including testing across multiple platforms and automatic publishing to PyPI.

## 🚀 Getting Started

Follow these steps to start your new project from the private Gitea instance.

### 1. Create Project from Gitea Template

Instead of cloning, you should create a new repository directly from this template on your private Gitea instance.

1. Navigate to the template repository page on Gitea.

2. Click the **"Use this template"** button.

3. Fill in the details for your new repository and create it.

4. Clone your newly created repository to your local machine.

```bash

# Clone the repository YOU created from the template

git clone <your-new-repository-url>

cd <your-new-repository>

```

### 2. Rename the Project

This is the most important first step. We've provided a script to automate this.

Run the following command and enter your new package name when prompted (e.g., `my_awesome_app`). Use underscores (`_`) for package names, not hyphens (`-`).

```bash

make rename

```

This command will:

- Rename the `src/myrepositorytemplate` directory to `src/your_new_name`.

- Update the `name` in `pyproject.toml`.

- Update the `PACKAGE_NAME` variable in the `Makefile`.

### 3. Initialize the Environment & Activate

Initialize your development environment with this interactive command.

```bash
make init
```

This powerful command configures your project by:

1. **Checking for `uv`**: If `uv` is not found, it will be automatically installed.
2. **Prompting for Python Version**: It will ask you to enter a Python version for the project (e.g., `3.11`, `3.12`).
3. **Pinning the Version**: It saves your choice to a `.python-version` file and updates the `requires-python` field in `pyproject.toml`.
4. **Installing Dependencies**: It creates a virtual environment and installs all dependencies from `uv.lock`.

Finally, activate the environment to begin working:

```bash
source .venv/bin/activate
```

You are now ready to code!

---

## 📖 Command Reference

All common tasks are managed through `make` commands. Here is a comprehensive reference.

### 📦 Dependency Management

Manage your project's dependencies with these commands.

- `make init` or `make install`: The primary **interactive** command to set up the project. It pins the Python version, updates metadata, and installs all dependencies.

- `make install-uv`: A helper command to install the `uv` package manager itself, which is called automatically by `make init` if needed.

- `make sync`: Use this if you have manually changed `pyproject.toml` or pulled new changes. It syncs the virtual environment with the `uv.lock` file.

- `make update`: Upgrades all dependencies to the latest allowed versions according to `pyproject.toml` and updates the `uv.lock` file.

- `make tree`: Displays the complete dependency tree, useful for debugging dependency conflicts.

To add or remove dependencies, use `uv` directly:

- `uv add <package>`: Add a new main dependency.

- `uv add --dev <package>`: Add a new development dependency.

- `uv remove <package>`: Remove a dependency.

### 🎨 Code Quality

Ensure your code stays clean, formatted, and type-safe.

- `make format`: Formats all code in the project using `Ruff Formatter`.

- `make lint`: Lints all code using `Ruff Linter` and attempts to auto-fix any fixable issues.

- `make type-check`: Runs `Mypy` to perform a strict static type analysis of your code.

- `make check`: The all-in-one quality command. It runs `format`, `lint`, and `type-check` sequentially. **Run this before every commit!**

### 🧪 Testing

Run your test suite and check for code coverage.

- `make test`: Executes the entire test suite using `pytest`.

- `make test-cov`: Runs the tests and generates a detailed HTML coverage report. It will then automatically open the report in your default web browser for inspection.

### 🏗️ Build & Release

This template is configured with a dual-release strategy for internal and public releases.

- `make build`: A safe build command. It first runs `make clean`, `make check`, and `make test`. If everything passes, it builds the distribution packages (`.whl` and `.tar.gz`) into the `dist/` directory.

- `make release`: This is a helper command that **does not perform a release**. Instead, it prints out the URLs and instructions for performing a manual release on Gitea or GitHub, reminding you of the process.

#### Internal Development Snapshots (Gitea)

- **Trigger**: `git push` to the `main` branch on the **private Gitea repository**.

- **Outcome**: Gitea Actions automatically builds a development snapshot (e.g., `0.1.0.dev5`) and publishes it to the **internal Gitea package registry**.

#### Public Stable Releases (GitHub)

- **Trigger**: Manually creating a new **Release** on the **public GitHub mirror**.

- **Outcome**: GitHub Actions runs extensive tests, builds a stable version (e.g., `0.1.0`), publishes it to **public PyPI**, and attaches the artifacts to the GitHub Release page.

### 🧹 Utilities & Maintenance

Keep your project directory clean.

- `make clean`: Removes temporary files and build artifacts, such as `__pycache__`, `.pytest_cache`, `build/`, `dist/`, and coverage reports.

- `make distclean`: A more thorough cleaning. It performs `make clean` and then **completely removes the virtual environment (`.venv`) and `.python-version` file**. Use this when you want to reset your environment from scratch.

### 🐳 Docker Integration

This template includes pre-configured Docker support for development and deployment.

- `make docker-build`: Builds the Docker image for your application. This is useful for testing the image build process or preparing for deployment.

- `make docker-up`: Builds (if necessary) and starts your application services in detached mode using `docker compose`. This is the primary way to run your application in a Dockerized environment during development.

- `make docker-down`: Stops and removes the containers, networks, and volumes created by `docker compose up`. Use this to clean up your Docker environment when you're done.

You can view the running container logs with `docker compose logs -f` after running `make docker-up`.
