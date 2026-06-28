# =============================================================================
# 🛠️ MyRepositoryTemplate Makefile (Production Ready)
# =============================================================================

# ⚠️ 模板使用者请修改这里：你的包名（对应 src/ 下的目录名）
PACKAGE_NAME := molgr

# --- 自动检测版本号 ---
# 尝试通过 importlib 读取已安装包的版本，如果失败则显示 "dynamic"
VERSION := $(shell uv run python -c "from importlib.metadata import version; print(version('$(PACKAGE_NAME)'))" 2>/dev/null || echo "dynamic")
EXTENSION_MODULES := molgr._core
CPP_BUILD_DIR ?= build
# 检测操作系统，用于打开浏览器命令
DETECTED_OS := $(shell uname)
ifeq ($(DETECTED_OS), Darwin)
	OPEN_CMD := open
else
	OPEN_CMD := xdg-open
endif

.PHONY: help init install-uv install sync update tree format lint type-check check test test-cov clean distclean build release rename docker-build docker-up docker-down cpp-dev-install cpp-build cpp-ide-init stubs clean-cpp

# =============================================================================
# 📝 帮助文档
# =============================================================================
help:
	@echo "📚 \033[1;34mPython Project Makefile Helper\033[0m"
	@echo ""
	@echo "📦 \033[1;33mDependency Management:\033[0m"
	@echo "  make install-uv  ⬇️ Install uv (The package manager)"
	@echo "  make init        🚀 Initialize environment (Interactive: Pin Python + Update Metadata)"
	@echo "  make update      🔄 Update dependencies (uv lock --upgrade)"
	@echo "  make tree        🌳 Show dependency tree"
	@echo ""
	@echo "🎨 \033[1;33mCode Quality:\033[0m"
	@echo "  make format      ✨ Format code (ruff)"
	@echo "  make lint        🔍 Lint code (ruff check --fix)"
	@echo "  make type-check  🦆 Static type check (mypy)"
	@echo "  make check       🛡️ Run all checks (format + lint + type-check)"
	@echo ""
	@echo "🧪 \033[1;33mTesting:\033[0m"
	@echo "  make test        🌡️ Run unit tests"
	@echo "  make test-cov    📊 Run tests with HTML coverage report & open it"
	@echo ""
	@echo "🏗️ \033[1;33mBuild & Release:\033[0m"
	@echo "  make build       📦 Build package (sdist + wheel)"
	@echo "  make release     🚀 Print release instructions (Git Tag / GitHub Release)"
	@echo ""
	@echo "🧹 \033[1;33mUtilities:\033[0m"
	@echo "  make clean       🧹 Clean build artifacts & cache"
	@echo "  make distclean   🗑️ Clean EVERYTHING (including .venv)"
	@echo "  make rename      🏷️ [Template Tool] Global rename of '$(PACKAGE_NAME)'"
	@echo ""
	@echo "🐳 \033[1;33mDocker Integration:\033[0m"
	@echo "  make docker-build  🏗️ Build Docker image"
	@echo "  make docker-up     🚀 Run with Docker Compose"
	@echo "  make docker-down   🛑 Stop Docker services"
	@echo ""
	@echo "🦕 \033[1;33mC++ Dev:\033[0m"
	@echo "  make cpp-dev-install  📦 Install C++ dev dependencies"
	@echo "  make cpp-build        🔨 Rebuild editable C++ extension and refresh IDE config"
	@echo "  make cpp-ide-init     🧠 Generate compile_commands / clangd / VSCode C++ IDE config"
	@echo "     Optional: make cpp-ide-init CPP_BUILD_DIR=cmake-build-debug"
	@echo ""
	@echo "📌 Current Package: $(PACKAGE_NAME)"
	@echo "📌 Detected Version: $(VERSION)"

# =============================================================================
# 📦 依赖管理
# =============================================================================

# --- 🛠️ 安装 uv 工具 ---
install-uv:
	@echo "⬇️ Installing uv via official script..."
	@curl -LsSf https://astral.sh/uv/install.sh | sh
	@echo "✅ uv installed! You might need to restart your shell or run 'source $$HOME/.cargo/env'"

# --- 🚀 初始化项目 ---
init:
	@# 1. 检查 uv 是否存在
	@if ! command -v uv >/dev/null 2>&1; then \
		echo "⚠️ uv not found. Installing now..."; \
		$(MAKE) install-uv; \
	fi
	
	@# 2. 交互式询问 Python 版本
	@echo "🐍 \033[1;33mLet's configure your Python environment.\033[0m"
	@read -p "👉 Enter Python version to use (default: 3.11): " py_ver; \
	if [ -z "$$py_ver" ]; then \
		py_ver="3.11"; \
	fi; \
	\
	echo "📌 Pinning Python version to $$py_ver (.python-version)..."; \
	uv python pin $$py_ver; \
	\
	echo "📝 Updating pyproject.toml (requires-python >= $$py_ver)..."; \
	sed -i.bak "s/^requires-python = .*/requires-python = \">=$$py_ver\"/" pyproject.toml && rm pyproject.toml.bak

	@# 3. 同步依赖
	@echo "🚀 Installing dependencies..."
	uv sync --all-extras --dev
	
	@echo "✅ Environment ready! Activate with: source .venv/bin/activate"

# 别名
install: init

# 仅仅同步
sync:
	uv sync

# 升级依赖
update:
	@echo "🔄 Updating dependencies..."
	uv lock --upgrade
	uv sync

# 显示依赖树
tree:
	uv tree

# =============================================================================
# 🎨 代码质量
# =============================================================================
format:
	@echo "🎨 Running Ruff Formatter..."
	uv run ruff format .

lint:
	@echo "🔍 Running Ruff Linter..."
	uv run ruff check . --fix

type-check:
	@echo "🦆 Running Mypy Type Checker..."
	uv run mypy src

check: format lint type-check

# =============================================================================
# 🧪 测试与覆盖率
# =============================================================================
test:
	@echo "🧪 Running Pytest..."
	uv run pytest

test-cov:
	@echo "📊 Running Test Coverage..."
	uv run pytest --cov=src --cov-report=html
	@echo "🌍 Opening coverage report..."
	@$(OPEN_CMD) htmlcov/index.html

# =============================================================================
# 🏗️ 构建与打包
# =============================================================================
build: clean check test
	@echo "🏗️ Building package (Hatchling + UV)..."
	uv build

# =============================================================================
# 🚀 发布流程
# =============================================================================
release:
	@echo ""
	@echo "🚀 \033[1;32mReady to release version?\033[0m (Current detection: $(VERSION))"
	@echo "---------------------------------------------------"
	@echo "Since you are using CI/CD driven releases:"
	@echo "1. Commit all your changes."
	@echo "2. Go to your repository release page:"
	@echo "   🔗 GitHub: https://github.com/YourUser/$(PACKAGE_NAME)/releases/new"
	@echo "   🔗 Gitea:  (Your Gitea URL)/$(PACKAGE_NAME)/releases/new"
	@echo "3. Draft a new release and create a tag (e.g., v0.1.0)."
	@echo "---------------------------------------------------"

# =============================================================================
# 🧹 清理
# =============================================================================
distclean: clean
	@echo "🗑️ Removing virtual environment (.venv)..."
	rm -rf .venv .python-version
	@echo "✨ Project is clean. Run 'make init' to restart."

# =============================================================================
# 🛠️ 模板实用工具
# =============================================================================
rename:
	@read -p "Enter new package name (e.g., my_awesome_tool): " new_name; \
	if [ -z "$$new_name" ]; then \
		echo "❌ Name cannot be empty"; \
		exit 1; \
	fi; \
	if [ -d "src/$$new_name" ]; then \
		echo "❌ Directory src/$$new_name already exists"; \
		exit 1; \
	fi; \
	\
	echo "🔄 Renaming directory src/$(PACKAGE_NAME) -> src/$$new_name..."; \
	mv src/$(PACKAGE_NAME) src/$$new_name; \
	\
	echo "🔄 Replacing '$(PACKAGE_NAME)' with '$$new_name' in ALL project files..."; \
	find . -type f \
		-not -path "./.git/*" \
		-not -path "./.venv/*" \
		-not -path "./.mypy_cache/*" \
		-not -path "./.ruff_cache/*" \
		-not -path "./.pytest_cache/*" \
		-not -path "./__pycache__/*" \
		-not -path "*/__pycache__/*" \
		-not -path "./dist/*" \
		-not -path "./build/*" \
		-not -path "./uv.lock" \
		-exec grep -Iq "$(PACKAGE_NAME)" {} \; -print0 | \
		xargs -0 sed -i.bak "s/$(PACKAGE_NAME)/$$new_name/g"; \
	\
	echo "🧹 Cleaning up backup files..."; \
	find . -name "*.bak" -type f -delete; \
	\
	echo "✅ Rename complete!"; \
	echo "👉 Note: 'uv.lock' might be out of sync. Please run 'make init' or 'uv sync' to regenerate it."

# =============================================================================
# 🐳 Docker 常用命令
# =============================================================================
docker-build:
	@echo "🏗️ Building Docker image for $(PACKAGE_NAME)..."
	docker build -t $(PACKAGE_NAME):latest .

docker-up:
	@echo "🚀 Starting services..."
	docker compose up -d --build
	@echo "📜 Use 'docker compose logs -f' to follow logs."

docker-down:
	@echo "🛑 Stopping services..."
	docker compose down

# =============================================================================
# 🦕 C++ 扩展开发 (🆕 新增板块)
# =============================================================================

cpp-dev-install:
	@echo "Installing C++ dev dependencies"
	uv pip install scikit-build-core pybind11 pybind11-stubgen setuptools-scm
	uv run python -c "import sys; raise SystemExit(0 if sys.version_info < (3, 10) else 1)" \
		&& uv pip install "openbabel-wheel>=3.1.1" \
		|| uv pip install "openbabel>=3.2.0"
# 🔨 快速编译 C++ 扩展
# 使用 --no-build-isolation 避免每次重新安装构建依赖，加快重编速度
# -v 显示 CMake/编译器 输出，方便调试
cpp-build:
	@echo "🔨 Re-compiling C++ extension (Editable Mode)..."
	uv pip install -e . -v --no-build-isolation
	@$(MAKE) cpp-ide-init CPP_BUILD_DIR=$(CPP_BUILD_DIR)

cpp-ide-init:
	@echo "🧠 Generating C++ IDE configuration..."
	uv run python scripts/gen_vscode_config_with_ob.py --build-dir $(CPP_BUILD_DIR)

# 🤖 生成类型存根 (Stubs)
# 1. 确保安装了 pybind11-stubgen
# 2. 运行生成器，输出到 src/ 目录以便编辑器识别
# --root-module-suffix="" 是为了避免生成 molgr-stubs 文件夹，而是直接在 src/molgr 下生成 .pyi
stubs:
	@echo "🤖 Checking for pybind11-stubgen..."
	@uv run python -c "import pybind11_stubgen" 2>/dev/null || (echo "⚠️ pybind11-stubgen not found. Installing..." && uv pip install pybind11-stubgen)
	@echo "🧹 Removing stale C++ backend stubs..."
	@rm -rf src/molgr/_core/dev src/molgr/_core/dev.pyi
	@rm -f src/molgr/_core/pipeline.pyi
	@rm -rf src/molgr/_core/stages src/molgr/_core/stages.pyi
	@echo "🤖 Generating Type Stubs for: $(EXTENSION_MODULES)..."
	@for module in $(EXTENSION_MODULES); do \
		echo "   Processing $$module..."; \
		uv run pybind11-stubgen $$module \
			--output-dir src/ \
			--root-suffix "" \
			--ignore-all-errors \
			--numpy-array-use-type-var; \
	done
	@if [ -f src/molgr/_core/pipeline.pyi ]; then \
		mkdir -p src/molgr/_core/pipeline; \
		mv src/molgr/_core/pipeline.pyi src/molgr/_core/pipeline/__init__.pyi; \
	fi
	@uv run python scripts/split_pipeline_stub.py
	@echo "🎨 Formatting generated stubs..."
	@uv run ruff format src/molgr/_core
	@echo "🔍 Linting generated stubs..."
	@uv run ruff check src/molgr/_core --fix
	@echo "✅ Stubs generated in src/!"

# =============================================================================
# 🧹 清理 (🔄 修改板块)
# =============================================================================
clean:
	@echo "🧹 Cleaning artifacts..."
	rm -rf dist build htmlcov coverage.xml .coverage
	rm -f compile_commands.json
	# 🆕 清理 scikit-build-core 的构建缓存
	rm -rf _skbuild
	# 🆕 清理编译出来的二进制扩展 (.so, .pyd)
	find src -type f -name "*.so" -delete
	find src -type f -name "*.pyd" -delete
	# 清理 Python 缓存
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name ".ruff_cache" -exec rm -rf {} +
	find . -type d -name ".mypy_cache" -exec rm -rf {} +
