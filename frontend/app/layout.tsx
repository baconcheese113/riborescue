import "@fontsource-variable/inter";
import "./globals.css";
import type { ReactNode } from "react";
import { Nav } from "./nav";

export const metadata = {
  title: "RiboRescue",
  description: "Matching nonsense variants to candidate readthrough therapies.",
};

export default function RootLayout({ children }: { children: ReactNode }) {
  return (
    <html lang="en">
      <body>
        <Nav />
        {children}
      </body>
    </html>
  );
}
