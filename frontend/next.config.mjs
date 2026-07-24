/** Static export, so the app is plain files on Cloudflare Pages and never runs a server. */
const nextConfig = {
  output: "export",
  images: { unoptimized: true },
};

export default nextConfig;
