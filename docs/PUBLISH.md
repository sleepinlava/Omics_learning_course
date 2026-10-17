# 发布到公开 GitHub 仓库

目标：`sleepinlava/ximcs`，公开。当前环境无法连接 GitHub API，仓库尚未创建；本地源文件已准备，可使用 `dist/ximcs.bundle` 恢复完整提交。

在可访问 GitHub、且 `gh auth status` 已确认登录为有权创建该仓库的账号的终端运行：

```bash
# 从本工作区根目录运行；使用全新的目标目录，避免覆盖已有文件。
git clone -b main dist/ximcs.bundle ../ximcs-publish
cd ../ximcs-publish
# clone bundle 会产生一个指向本地 bundle 的 origin，先为它改名。
git remote rename origin local-bundle
gh repo create sleepinlava/ximcs --public --source=. --remote=origin --push
```

若远程仓库已由其他步骤创建，先核对其归属和内容，不强制覆盖。空仓库可改用 `git remote add origin https://github.com/sleepinlava/ximcs.git` 与 `git push -u origin main`。非空仓库先 fetch/review 再整合。

`dist/ximcs-source.tar.gz` 是同一提交的源码归档，`dist/SHA256SUMS` 用于校验。不要上传本地 `data/raw_fastq/`、第三方完整教程、依赖环境或所有历史输出。成功后用 `gh repo view sleepinlava/ximcs --json url,visibility` 和 `git ls-remote origin refs/heads/main` 核对公开状态与提交号。
