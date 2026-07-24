"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

// The four routes share one bar. Active state is drawn from the path, links are real anchors so the
// keyboard and the browser's history both work, and the whole thing wraps on a narrow screen.
const ROUTES: [string, string][] = [
  ["/", "Home"],
  ["/patient", "Patient"],
  ["/researcher", "Researcher"],
  ["/explorer", "Explorer"],
];

export function Nav() {
  const pathname = usePathname();
  return (
    <nav className="nav" aria-label="Primary">
      <Link href="/" className="brand">
        RiboRescue
      </Link>
      <div className="nav-links">
        {ROUTES.slice(1).map(([href, label]) => {
          const active = pathname === href;
          return (
            <Link
              key={href}
              href={href}
              className={active ? "nav-link active" : "nav-link"}
              aria-current={active ? "page" : undefined}
            >
              {label}
            </Link>
          );
        })}
      </div>
    </nav>
  );
}
