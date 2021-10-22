import throttle from "lodash.throttle";
import * as React from "react";
import * as ReactDOM from "react-dom";
import {CollectionContext} from '../CollectionContext';

export function useSyncedRef<T>(value: T): React.MutableRefObject<T> {
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


export function useNameInput(id: string, initialValue: string): [string, React.Dispatch<React.SetStateAction<string>>, React.ReactNode] {
  const [name, setName] = React.useState<string>("");
  const uniqueRef = React.useRef(id);
  const collections = React.useContext(CollectionContext);

  // Remove any [ ] . characters as they are matched with lodash.get and break lineup for example.
  const currentName = (name || initialValue).replaceAll(/\[|\]|\./g, '_');

  const nameAlreadyTaken = collections.some((c) => c.name === currentName);

  return [currentName, setName, <div className="mb-3">
    <label htmlFor={uniqueRef.current}>Name</label>
    <input
      type="text"
      className={`form-control form-control-sm ${nameAlreadyTaken ? "is-invalid" : ""}`}
      id={uniqueRef.current}
      required={nameAlreadyTaken ? true : !initialValue} // Required depends if we have a valid placeholder
      pattern={nameAlreadyTaken ? '^$' : undefined} // Add a regex matching only the empty string if the name is already taken
      value={name}
      placeholder={initialValue}
      onChange={(e) => setName(e.currentTarget.value)}
    />
    {nameAlreadyTaken ? (
      <div className="invalid-feedback">Collection with this name already exists.</div>
    ) : null}
  </div>];
}

export function timeIt<T>(name: string, f: () => T): T {
  console.time(name);
  const result = f();
  console.timeEnd(name);
  return result;
}