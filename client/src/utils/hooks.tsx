import throttle from "lodash.throttle";
import * as React from "react";
import * as ReactDOM from "react-dom";

export function useSyncedRef<T>(value: T) {
  const valueRef = React.useRef<T>(value);
  valueRef.current = value;
  return valueRef;
}

export const useMousePosition = (active: boolean): { x: number; y: number } | null => {
  const [position, setPosition] = React.useState<{
    x: number;
    y: number;
  } | null>(null);

  React.useEffect(() => {
    if (active) {
      const setFromEvent = throttle((e: MouseEvent) => setPosition({ x: e.clientX, y: e.clientY }), 5);
      window.addEventListener("mousemove", setFromEvent);

      return () => {
        window.removeEventListener("mousemove", setFromEvent);
      };
    } else {
      setPosition(null);
    }
  }, [active]);

  return position;
};

export const Tooltip = ({ children, anchor = 'right' }: { children: React.ReactNode, anchor?: 'top' | 'right' }) => {
  const [node, setNode] = React.useState<HTMLDivElement | null>(null);
  const position = useMousePosition(Boolean(node && children));

  React.useEffect(() => {
    const container = document.createElement("div");
    container.style.position = "absolute";
    container.style.pointerEvents = "none";
    container.style.transform = anchor === 'right' ? "translate(10%, -50%)" : "translate(-50%, 5%)";
    container.style.zIndex = "10000";
    document.body.appendChild(container);
    setNode(container);

    return () => {
      document.body.removeChild(container);
    };
  }, [anchor]);

  React.useEffect(() => {
    if (node) {
      if (position) {
        node.style.left = `${position.x}px`;
        node.style.top = `${position.y}px`;
        node.style.display = "block";
      } else {
        node.style.display = "none";
      }
    }
  }, [node, position]);

  return node && children && position ? ReactDOM.createPortal(children, node) : null;
};

export function usePrevious<T>(value: T | null): T | null {
  // The ref object is a generic container whose current property is mutable ...
  // ... and can hold any value, similar to an instance property on a class
  const ref = React.useRef<T | null>(null);

  // Store current value in ref
  React.useEffect(() => {
    ref.current = value;
  }, [value]); // Only re-run if value changes

  // Return previous value (happens before update in useEffect above)
  return ref.current;
}
