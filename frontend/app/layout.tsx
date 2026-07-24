import "./globals.css";
import type { ReactNode } from "react";

export const metadata = {
  title: "RiboRescue — variant × therapy",
  description: "Which readthrough therapies a nonsense variant might be amenable to.",
};

export default function RootLayout({ children }: { children: ReactNode }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
