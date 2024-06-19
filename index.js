const { execSync } = require("child_process");
const fs = require("fs");
const path = require("path");

const files = [];

const listFiles = (dir) => {
  fs.readdirSync(dir).forEach((file) => {
    // Do not recurse on node_modules
    if (file === "node_modules" || file === ".git") return;

    // Add file to list
    const filePath = path.join(dir, file);
    const stat = fs.statSync(filePath);
    if (stat.isDirectory()) {
      listFiles(filePath);
    } else {
      files.push({
        filePath: path.relative(__dirname, filePath).replace(/[/\\]/, "/"),
        lastModified: stat.mtime,
      });
    }
  });
};

listFiles(__dirname);

files.sort((a, b) => (a.lastModified < b.lastModified ? -1 : 1));

const datestring = (d) =>
  ("0" + d.getDate()).slice(-2) +
  "-" +
  ("0" + (d.getMonth() + 1)).slice(-2) +
  "-" +
  d.getFullYear() +
  " " +
  ("0" + d.getHours()).slice(-2) +
  ":" +
  ("0" + d.getMinutes()).slice(-2) +
  ":" +
  ("0" + d.getSeconds()).slice(-2);

const createCommit = async () => {
  for (const f of files) {
    if (f.filePath.includes(" ")) {
      const newFilePath = f.filePath.replace(/ /g, "_");
      fs.renameSync(
        path.resolve(__dirname, f.filePath),
        path.resolve(__dirname, newFilePath)
      );
      f.filePath = newFilePath;
    }

    execSync(`git add "${f.filePath}"`);
    await new Promise((resolve) => setTimeout(resolve, 100));
    execSync(
      `git commit --date "${datestring(f.lastModified)}" -m "Create ${
        f.filePath
      }"`
    );
  }
};

createCommit();
