# GitHub Pages Setup Guide for pyWC

Your MkDocs documentation site has been created and committed. Follow these steps to enable GitHub Pages:

## Step 1: Enable GitHub Pages

1. Go to your repository: https://github.com/octupole/pyWC

2. Click on **Settings** (top tab)

3. In the left sidebar, click **Pages**

4. Under "Build and deployment":
   - **Source**: Select **Deploy from a branch**
   - **Branch**: Select **gh-pages** (this will be created by the GitHub Action)
   - **Folder**: Leave as **/ (root)**

5. Click **Save**

## Step 2: Wait for GitHub Actions

The documentation will be automatically built and deployed when you push to `main`.

Check the progress:
1. Go to the **Actions** tab in your repository
2. Look for the "Deploy Documentation" workflow
3. Wait for it to complete (green checkmark)

## Step 3: Access Your Documentation Site

Once the workflow completes, your documentation will be live at:

**https://octupole.github.io/pyWC/**

## Step 4: Add Website Link to Repository

1. Go to your repository homepage
2. Click the **⚙️ gear icon** next to "About" (top right)
3. Check **"Use your GitHub Pages website"**
4. Click **Save changes**

This will add a link to your documentation in the repository sidebar.

## Troubleshooting

### GitHub Action Fails

If the "Deploy Documentation" action fails:

1. Check the Actions tab for error messages
2. Common issues:
   - **Permissions**: Go to Settings > Actions > General > Workflow permissions
     - Select **"Read and write permissions"**
     - Check **"Allow GitHub Actions to create and approve pull requests"**
   - **Branch doesn't exist**: The first run creates `gh-pages` branch automatically

### Documentation Not Updating

1. Check that you pushed to the `main` branch
2. Verify the GitHub Action ran successfully
3. Clear your browser cache (Ctrl+F5 or Cmd+Shift+R)
4. Wait a few minutes for GitHub's CDN to update

### 404 Error

If you get a 404:
1. Ensure the `gh-pages` branch was created
2. Check Settings > Pages shows the correct source
3. Wait 5-10 minutes after first deployment

## Building Documentation Locally

To preview the documentation before pushing:

```bash
# Install dependencies
pip install -r docs/requirements.txt

# Serve locally
mkdocs serve
```

Then visit http://127.0.0.1:8000

## Updating Documentation

1. Edit files in `docs/` directory
2. Preview with `mkdocs serve`
3. Commit and push to `main`
4. GitHub Actions automatically deploys

## Custom Domain (Optional)

To use a custom domain like `docs.yoursite.com`:

1. Add a CNAME record in your DNS:
   ```
   docs.yoursite.com  →  octupole.github.io
   ```

2. In GitHub Settings > Pages:
   - Enter your custom domain in "Custom domain"
   - Click Save
   - Wait for DNS check to complete

3. Create `docs/CNAME` file:
   ```
   docs.yoursite.com
   ```

4. Commit and push

## Next Steps

Once GitHub Pages is enabled:

1. ✅ Visit your documentation site
2. ✅ Share the link with users
3. ✅ Add badges to README.md:
   ```markdown
   [![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://octupole.github.io/pyWC/)
   ```
4. ✅ Tweet/announce your documentation

## Support

If you encounter issues:
- GitHub Pages documentation: https://docs.github.com/en/pages
- MkDocs Material: https://squidfunk.github.io/mkdocs-material/
- Open an issue: https://github.com/octupole/pyWC/issues

---

**Your documentation is ready to go! Just enable GitHub Pages in the repository settings.**
