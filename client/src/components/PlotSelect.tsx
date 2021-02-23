import * as React from "react";

interface GenericOption<T> {
  option: T | null | undefined;
  setOption(option: T | null): void;
}

export function PlotSelect({
  multi,
  option,
  setOption,
  label,
  options,
  id,
  disabled,
}: {
  options: string[];
  disabled?: boolean;
  id: string;
  label?: string;
} & (
  | ({
      multi: true;
    } & GenericOption<string[]>)
  | ({
      multi?: false;
    } & GenericOption<string>)
)) {
  // @ts-ignore
  const values: (string | null)[] = (multi ? option || [] : [option]).filter(Boolean);
  if (values.length === 0 || (multi && values[values.length - 1] !== null)) {
    values.push(null);
  }

  return (
    <div className="d-flex" style={{ opacity: values.filter(Boolean).length === 0 ? 0.8 : undefined }}>
      {label ? (
        <label className="my-1 mr-2" style={{ whiteSpace: "nowrap" }} htmlFor={id}>
          {label}
        </label>
      ) : null}
      <div className="input-group my-1 mr-sm-2">
        {values.map((v, i) => (
          <select
            key={i}
            id={id}
            className="custom-select custom-select-sm"
            disabled={options.length === 0 || disabled}
            value={v || ""}
            onChange={(e) => {
              const currentValue = e.currentTarget.value || null;
              if (!multi) {
                setOption(currentValue as any);
              } else {
                const newValues = currentValue ? [...values] : values.slice(0, i);
                if (currentValue) {
                  newValues[i] = currentValue;
                }
                setOption(newValues.filter(Boolean) as any);
              }
            }}
          >
            <option value={""}>{options.length === 0 ? "No available properties" : "Choose..."}</option>
            {options.map((key) => (
              <option key={key} value={key}>
                {key}
              </option>
            ))}
          </select>
        ))}
      </div>
    </div>
  );
}
