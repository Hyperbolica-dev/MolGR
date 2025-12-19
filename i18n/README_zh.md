# 开发者指南: Modern Python Project Template

[Read this guide in English](../README_DEV.md)

本文档面向使用此项目模板的开发者，旨在说明构建高质量 Python 应用的理念、工具和工作流。

## 核心技术栈

本模板构建于一系列现代、高效的工具基础之上：

- **`uv`**: 一个速度极快的 Python 包安装器和解析器，用于所有依赖管理、环境创建和任务执行。
- **`Ruff`**: 一个速度极快的 Python Linter 和代码格式化工具。
- **`Mypy`**: 事实上的静态类型检查标准。
- **`Pytest`**: 用于编写健壮、可扩展测试的框架。
- **`Hatch`**: 作为项目构建后端，版本控制由 `hatch-vcs` 管理。
- **GitHub Actions**: 用于自动化 CI/CD，包括跨平台测试和向 PyPI 自动发布。

## 🚀 快速上手

请遵循以下步骤在私有 Gitea 实例上启动您的新项目。

### 1. 从 Gitea 模板创建项目

您应该直接在私有的 Gitea 实例上从此模板创建一个新仓库，而不是克隆本项目。

1. 在 Gitea 上打开此模板仓库的页面。
2. 点击 **"使用此模板"** (Use this template) 按钮。
3. 填写您新仓库的详细信息并创建它。
4. 将您 **新创建的仓库** 克隆到本地。

```bash
# 克隆您从模板创建的新仓库
git clone <your-new-repository-url>
cd <your-new-repository>
```

### 2. 重命名项目

这是最关键的第一步。我们提供了一个自动化脚本来完成此操作。

运行以下命令，并在提示时输入您的新包名（例如 `my_awesome_app`）。包名请使用下划线（`_`），而不是连字符（`-`）。

```bash
make rename
```

此命令将：

- 将 `src/myrepositorytemplate` 目录重命名为 `src/your_new_name`。
- 更新 `pyproject.toml` 中的 `name` 字段。
- 更新 `Makefile` 中的 `PACKAGE_NAME` 变量。

### 3. 初始化并激活环境

使用这个交互式命令来初始化您的开发环境。

```bash
make init
```

这个强大的命令会通过以下步骤配置您的项目：

1. **检查 `uv`**：如果您系统中没有 `uv`，它将被自动安装。
2. **提示输入 Python 版本**：它会请求您为项目输入一个 Python 版本（例如 `3.11`, `3.12`）。
3. **锁定版本**：它会将您的选择保存到 `.python-version` 文件，并更新 `pyproject.toml` 中的 `requires-python` 字段。
4. **安装依赖**：它会创建一个虚拟环境，并从 `uv.lock` 文件中安装所有依赖。

最后，激活环境以开始工作：

```bash
source .venv/bin/activate
```

现在您可以开始编码了！

---

## 📖 命令参考手册

所有常见任务都通过 `make` 命令管理。这是一份完整的参考。

### 📦 依赖管理

使用这些命令来管理您项目的依赖关系。

- `make init` 或 `make install`: 设置项目的主要 **交互式** 命令。它会锁定 Python 版本、更新元数据并安装所有依赖。
- `make install-uv`: 一个辅助命令，用于安装 `uv` 包管理器。如果需要，`make init` 会自动调用它。
- `make sync`: 当您手动修改了 `pyproject.toml` 或拉取了新变更时使用。它会根据 `uv.lock` 文件同步虚拟环境。
- `make update`: 根据 `pyproject.toml` 的规则，将所有依赖项升级到允许的最新版本，并更新 `uv.lock` 文件。
- `make tree`: 显示完整的依赖关系树，用于调试依赖冲突。

要添加或移除依赖，请直接使用 `uv`：

- `uv add <package>`: 添加一个新的主依赖。
- `uv add --dev <package>`: 添加一个新的开发依赖。
- `uv remove <package>`: 移除一个依赖。

### 🎨 代码质量

确保您的代码保持整洁、格式统一且类型安全。

- `make format`: 使用 `Ruff Formatter` 格式化项目中的所有代码。
- `make lint`: 使用 `Ruff Linter` 检查所有代码，并尝试自动修复可修复的问题。
- `make type-check`: 运行 `Mypy` 对您的代码执行严格的静态类型分析。
- `make check`: 一体化的质量检查命令。它会依次运行 `format`、`lint` 和 `type-check`。**请在每次提交前运行此命令！**

### 🧪 测试

运行您的测试套件并检查代码覆盖率。

- `make test`: 使用 `pytest` 执行完整的测试套件。
- `make test-cov`: 运行测试并生成一份详细的 HTML 格式的覆盖率报告。随后，它会自动在您的默认浏览器中打开该报告供您查阅。

### 🏗️ 构建与发布

本模板配置了用于内部和公开发布的双重发布策略。

- `make build`: 一个安全的构建命令。它会首先运行 `make clean`、`make check` 和 `make test`。如果一切顺利，它会将项目构建为分发包（`.whl` 和 `.tar.gz`）并存放在 `dist/` 目录中。
- `make release`: 这是一个辅助命令，它 **本身不执行发布**。相反，它会打印出在 Gitea 或 GitHub 上执行手动发布所需的 URL 和操作指南，以提醒您发布流程。

#### 内部开发快照 (Gitea)

- **触发条件**: 向 **私有 Gitea 仓库** 的 `main` 分支执行 `git push`。
- **最终结果**: Gitea Actions 自动构建一个开发快照版本 (例如 `0.1.0.dev5`) 并将其发布到 **内部 Gitea 包注册中心**。

#### 公开稳定版本 (GitHub)

- **触发条件**: 在 **公共 GitHub 镜像** 上手动创建一个新的 **Release**。
- **最终结果**: GitHub Actions 运行全面的测试，构建一个稳定版本 (例如 `0.1.0`)，将其发布到 **公共 PyPI**，并将软件包产物附加到 GitHub Release 页面。

### 🧹 工具与维护

保持您的项目目录干净整洁。

- `make clean`: 移除构建产物和临时文件，例如 `__pycache__`、`.pytest_cache`、`build/`、`dist/` 和覆盖率报告。
- `make distclean`: 更彻底的清理。它会先执行 `make clean`，然后 **彻底删除虚拟环境目录 (`.venv`) 和 `.python-version` 文件**。当您想从头开始重置您的环境时使用此命令。

### 🐳 Docker 集成

此模板包含预配置的 Docker 支持，用于开发和部署。

- `make docker-build`：构建你的应用程序的 Docker 镜像。这对于测试镜像构建过程或准备部署非常有用。

- `make docker-up`：使用 `docker compose` 构建（如果需要）并以分离模式启动你的应用程序服务。这是在开发过程中以 Docker 化环境运行应用程序的主要方式。

- `make docker-down`：停止并移除由 `docker compose up` 创建的容器、网络和卷。完成工作后使用此命令清理你的 Docker 环境。

运行 `make docker-up` 后，你可以使用 `docker compose logs -f` 查看正在运行的容器日志。
